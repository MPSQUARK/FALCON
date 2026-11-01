using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.Diagnostics;
using Falcon.Core.Services.IO;
using Falcon.Core.Services.Statistics;
using Serilog;

namespace Falcon.Core.Pipeline;

internal sealed class FireflyFitOrchestrator
{
    private readonly IDustService _dust;
    private readonly IStatisticsService _statistics;
    private readonly ISpectralFitter _fitter;
    private readonly IFireflyFitsWriter _fitsWriter;

    public FireflyFitOrchestrator(
        IDustService dust,
        IStatisticsService statistics,
        ISpectralFitter fitter,
        IFireflyFitsWriter fitsWriter)
    {
        _dust = dust;
        _statistics = statistics;
        _fitter = fitter;
        _fitsWriter = fitsWriter;
    }

    public bool TryFit(FireflyFitContext ctx, out FireflyFitOutput? output)
    {
        output = null;
        try
        {
            int nModels = ctx.ModelAges.Length;
            var modelWave = ctx.ModelWavelength;
            var modelFluxRaw = ctx.ModelFluxRaw;
            var dataWave = ctx.DataWavelength;
            var dataFlux = ctx.DataFlux;
            var dataError = ctx.DataError;
            var flags = Enumerable.Repeat(true, dataWave.Length).ToArray();

            double[] wave;
            double[] matchedData;
            double[] matchedError;
            double[][] matchedModels;

            if (ctx.ModelsOnDataGrid)
            {
                wave = dataWave;
                matchedData = dataFlux;
                matchedError = dataError;
                matchedModels = modelFluxRaw;
            }
            else
            {
                (wave, matchedData, matchedError, matchedModels) = _statistics.MatchDataModels(
                    dataWave, dataFlux, flags, dataError, modelWave, modelFluxRaw,
                    LegacyConfig.WaveLimits[0], LegacyConfig.WaveLimits[1]);
            }

            int nWave = Math.Min(wave.Length, Math.Min(matchedData.Length, Math.Min(matchedError.Length, matchedModels[0].Length)));
            wave = wave.Take(nWave).ToArray();
            matchedData = matchedData.Take(nWave).ToArray();
            matchedError = matchedError.Take(nWave).ToArray();
            for (int m = 0; m < matchedModels.Length; m++)
                matchedModels[m] = matchedModels[m].Take(nWave).ToArray();

            FitTrace.Emit("match", new Dictionary<string, object>
            {
                ["source"] = Path.GetFileName(ctx.DataFilePath),
                ["n_wave"] = nWave,
                ["wave_min"] = wave.Length > 0 ? wave[0] : 0,
                ["wave_max"] = wave.Length > 0 ? wave[^1] : 0,
                ["models_on_data_grid"] = ctx.ModelsOnDataGrid,
            }.Concat(FitTrace.Stats("data", matchedData))
             .Concat(FitTrace.Stats("error", matchedError))
             .Concat(FitTrace.ModelGridStats(matchedModels))
             .ToDictionary(kv => kv.Key, kv => kv.Value));

            double[][] normModels;
            double[] massFactors;
            using (FitTrace.BeginSpan("normalise"))
            {
                var (normModelsResult, massFactorsResult) = _statistics.NormaliseSpec(matchedData, matchedModels);
                normModels = normModelsResult;
                massFactors = massFactorsResult;
            }
            double bestEbv = 0;
            double[][] fitModels = normModels;
            double[] fitData = matchedData;
            double[] fitError = matchedError;

            int nanModels = matchedModels.Sum(m => m.Count(v => !double.IsFinite(v)));
            if (nanModels > 0)
                Log.Warning("Matched models contain {NanCount} non-finite values for {DataFile}", nanModels, ctx.DataFilePath);

            FitTrace.Emit("normalise", new Dictionary<string, object>
            {
                ["mass_factor_0"] = massFactors.Length > 0 ? massFactors[0] : 0,
                ["mass_factor_median"] = massFactors.Length > 0 ? Median(massFactors) : 0,
            }.Concat(FitTrace.ModelGridStats(normModels, "norm0"))
             .ToDictionary(kv => kv.Key, kv => kv.Value));

            if (LegacyConfig.HPF_Mode)
            {
                var dustInput = new DustAttenuationInput
                {
                    Wavelength = wave,
                    DataFlux = matchedData,
                    ErrorFlux = matchedError,
                    ModelFlux = matchedModels,
                    Age = ctx.ModelAges,
                    Metal = ctx.ModelMetals,
                    MaxEbv = LegacyConfig.MaxEbv,
                    NumDustVals = LegacyConfig.NumDustVals,
                    SmoothingLength = LegacyConfig.DustSmoothingLength,
                    DustLaw = DustLawName()
                };
                var dustResult = _dust.DetermineAttenuation(dustInput, _fitter);
                bestEbv = dustResult.Ebv;
                var attenuated = new double[nModels][];
                for (int m = 0; m < nModels; m++)
                {
                    attenuated[m] = new double[wave.Length];
                    for (int i = 0; i < wave.Length; i++)
                        attenuated[m][i] = matchedModels[m][i] * dustResult.AttenuationCurve[i];
                }
                (fitModels, massFactors) = _statistics.NormaliseSpec(matchedData, attenuated);
                FitTrace.Emit("dust", new Dictionary<string, object>
                {
                    ["ebv"] = bestEbv,
                    ["atten_median"] = Median(dustResult.AttenuationCurve),
                }.Concat(FitTrace.ModelGridStats(fitModels, "fit0"))
                 .ToDictionary(kv => kv.Key, kv => kv.Value));
            }
            else if (LegacyConfig.HpfOnlyMode)
            {
                var hpfData = _dust.Hpf(matchedData);
                var hpfModels = new double[nModels][];
                for (int m = 0; m < nModels; m++) hpfModels[m] = _dust.Hpf(matchedModels[m]);
                var zeroDat = new List<int>();
                for (int i = 0; i < hpfData.Length; i++)
                    if (double.IsNaN(hpfData[i]) || double.IsInfinity(hpfData[i]))
                    {
                        hpfData[i] = 0;
                        zeroDat.Add(i);
                    }
                for (int m = 0; m < nModels; m++)
                    foreach (var i in zeroDat) hpfModels[m][i] = 0;
                double medErr = Median(matchedError);
                double medData = Median(matchedData);
                double medHpf = Median(hpfData.Where(v => !double.IsNaN(v) && !double.IsInfinity(v)).DefaultIfEmpty(0).ToArray());
                fitError = new double[wave.Length];
                for (int i = 0; i < fitError.Length; i++)
                    fitError[i] = medData == 0 ? medErr : medErr / medData * medHpf;
                foreach (var i in zeroDat) fitError[i] = fitError.Max() * 999999.9;
                (fitModels, massFactors) = _statistics.NormaliseSpec(hpfData, hpfModels);
                fitData = hpfData;
            }

            var flatModels = FlattenModels(fitModels);
            var fitInput = new SpectralFitInput
            {
                Wavelength = wave,
                Data = fitData,
                Error = fitError,
                Models = flatModels
            };
            FitResult fitResult;
            using (FitTrace.BeginSpan("main_fit"))
            {
                fitResult = _fitter.Fit(fitInput, new FitterOptions
                {
                    MaxIterations = LegacyConfig.MaxIterations,
                    FitPerIterationCap = LegacyConfig.FitPerIterationCap,
                    MaxSspComponents = LegacyConfig.MaxSspComponents,
                    UseGpuEvaluator = LegacyConfig.UseGpuEvaluator,
                    ExperimentalMixSearch = LegacyConfig.ExperimentalMixSearch,
                    ExperimentalMixSearchMaxEvals = LegacyConfig.ExperimentalMixSearchMaxEvals,
                });
            }

            if (fitResult.Solutions.Count == 0)
            {
                Log.Warning("Firefly fit produced no solutions for {DataFile}", ctx.DataFilePath);
                return false;
            }

            var solutions = fitResult.Solutions
                .OrderBy(s => s.ChiSquared <= 0 ? double.MaxValue : s.ChiSquared)
                .ToList();
            var chis = solutions.Select(s => s.ChiSquared).ToArray();
            int dof = wave.Length;
            var probs = _statistics.ConvertChisToProbs(chis, dof);
            var lightWeights = solutions.Select(s => s.Weights).ToArray();
            var (unnormMass, massWeights) = _statistics.LightWeightsToMass(lightWeights, massFactors);

            if (massWeights.Any(row => row.Any(double.IsNaN)))
            {
                Log.Warning("Firefly fit produced NaN mass weights for {DataFile}", ctx.DataFilePath);
                return false;
            }

            double fluxUnits = Math.Pow(10, LegacyConfig.Flux_Scale_Factor);
            var averages = _statistics.CalculateAveragesPdf(
                probs, lightWeights, massWeights, unnormMass,
                ctx.ModelAges, ctx.ModelMetals, fluxUnits, ctx.DistanceLumCm, LegacyConfig.PdfSampling);

            var best = solutions[0];
            var bestFit = Blend(best.Weights, fitModels);
            double chiDofBest = dof > 0 && best.ChiSquared > 0 ? best.ChiSquared / dof : 0;

            FitTrace.Emit("fit", new Dictionary<string, object>
            {
                ["min_chi"] = fitResult.MinChiSquared,
                ["chi_dof_best"] = chiDofBest,
                ["n_solutions"] = solutions.Count,
                ["converged"] = fitResult.Converged,
                ["weight_0"] = best.Weights.Length > 0 ? best.Weights[0] : 0,
                ["weight_sum"] = best.Weights.Sum(),
                ["dof"] = dof,
            }.Concat(FitTrace.Stats("blend", bestFit))
             .ToDictionary(kv => kv.Key, kv => kv.Value));

            averages.TryGetValue("light_age", out double lightAge);
            averages.TryGetValue("light_metal", out double lightMetal);
            averages.TryGetValue("stellar_mass", out double stellarMass);
            FitTrace.Emit("pdf", new Dictionary<string, object>
            {
                ["light_age"] = lightAge,
                ["light_metal"] = lightMetal,
                ["stellar_mass"] = stellarMass,
                ["age_min_model"] = ctx.ModelAges.Length > 0 ? ctx.ModelAges.Min() : 0,
                ["age_max_model"] = ctx.ModelAges.Length > 0 ? ctx.ModelAges.Max() : 0,
            });

            var topSolutions = BuildTopSolutions(solutions, probs, fitModels, ctx.ModelAges, ctx.ModelMetals, dof, bestFit);

            string outDir = Path.Combine(LegacyConfig.PathOfProgram, "output");
            string baseName = Path.GetFileNameWithoutExtension(ctx.DataFilePath);
            string resultsPath = Path.Combine(outDir, baseName + ".results.json");
            string fitsPath = Path.Combine(outDir, baseName + "_firefly.fits");

            output = new FireflyFitOutput
            {
                SourceFile = ctx.DataFilePath,
                Converged = fitResult.Converged,
                BestEbv = bestEbv,
                MinChiSquared = fitResult.MinChiSquared,
                BestLightWeights = best.Weights,
                BestMassWeights = massWeights[0],
                Averages = averages,
                MatchedWave = wave,
                MatchedData = matchedData,
                MatchedError = matchedError,
                BestFitSpectrum = bestFit,
                OutputPath = fitsPath,
                ResultsPath = resultsPath,
                TopSolutions = topSolutions
            };

            if (LegacyConfig.Write_Results)
            {
                FalconResultsWriter.Write(output, ctx.Redshift);
                try
                {
                    _fitsWriter.WriteResults(fitsPath, output, ctx.Redshift);
                }
                catch (Exception ex)
                {
                    Log.Warning(ex, "FITS header sidecar write failed for {DataFile}; results JSON was written", ctx.DataFilePath);
                }
            }

            return true;
        }
        catch (Exception ex)
        {
            Log.Error(ex, "Firefly fit failed for {DataFile}", ctx.DataFilePath);
            return false;
        }
    }

    private static List<FalconSolutionRecord> BuildTopSolutions(
        List<FitSolution> orderedSolutions,
        double[] probs,
        double[][] fitModels,
        double[] agesGyr,
        double[] metalsLinear,
        int dof,
        double[]? bestFitSpectrum = null)
    {
        var records = new List<FalconSolutionRecord>();
        int rank = 0;
        for (int i = 0; i < orderedSolutions.Count && records.Count < LegacyConfig.TopResultsCount; i++)
        {
            var sol = orderedSolutions[i];
            rank++;
            double ageGyr = 0;
            double metalLin = 0;
            var components = new List<FalconSspComponent>();
            for (int m = 0; m < sol.Weights.Length; m++)
            {
                if (sol.Weights[m] <= 1e-8) continue;
                double ageG = agesGyr[m];
                ageGyr += sol.Weights[m] * ageG;
                metalLin += sol.Weights[m] * metalsLinear[m];
                components.Add(new FalconSspComponent
                {
                    ModelIndex = m,
                    LightWeight = sol.Weights[m],
                    AgeGyr = ageG,
                    Metal = metalsLinear[m] > 0 ? Math.Log10(metalsLinear[m]) : double.NaN
                });
            }

            double chi = sol.ChiSquared;
            if (!double.IsFinite(chi)) chi = 0;

            records.Add(new FalconSolutionRecord
            {
                Rank = rank,
                ChiSquared = chi,
                ChiSquaredPerDof = dof > 0 && double.IsFinite(chi) ? chi / dof : chi,
                Probability = i < probs.Length && double.IsFinite(probs[i]) ? probs[i] : 0,
                LightWeightedAgeGyr = ageGyr,
                LightWeightedMetal = metalLin,
                ModelFlux = rank == 1 && bestFitSpectrum is not null
                    ? bestFitSpectrum
                    : Blend(sol.Weights, fitModels),
                Components = components.OrderByDescending(c => c.LightWeight).ToList()
            });
        }
        return records;
    }

    private static string DustLawName() => LegacyConfig.Dust_Law switch
    {
        1 => "allen",
        2 => "prevot",
        _ => "calzetti"
    };

    private static double[] FlattenModels(double[][] models)
    {
        int n = models.Length;
        int w = models[0].Length;
        var flat = new double[n * w];
        for (int m = 0; m < n; m++)
            Array.Copy(models[m], 0, flat, m * w, w);
        return flat;
    }

    private static double[] Blend(double[] weights, double[][] models)
    {
        int len = models[0].Length;
        var result = new double[len];
        for (int m = 0; m < weights.Length; m++)
        {
            if (weights[m] == 0) continue;
            for (int i = 0; i < len; i++)
            {
                double v = models[m][i];
                if (double.IsFinite(v))
                    result[i] += weights[m] * v;
            }
        }
        return result;
    }

    private static double Median(double[] data)
    {
        var sorted = data.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n == 0 ? 0 : n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2 : sorted[n / 2];
    }
}

public sealed class FireflyFitContext
{
    public required string DataFilePath { get; init; }
    public required double[] DataWavelength { get; init; }
    public required double[] DataFlux { get; init; }
    public required double[] DataError { get; init; }
    public required double[] ModelWavelength { get; init; }
    public required double[][] ModelFluxRaw { get; init; }
    /// <summary>When true, model flux is already sampled on the data wavelength grid.</summary>
    public bool ModelsOnDataGrid { get; init; }
    public required double[] ModelAges { get; init; }
    public required double[] ModelMetals { get; init; }
    public required double DistanceLumCm { get; init; }
    public required float Redshift { get; init; }
}
