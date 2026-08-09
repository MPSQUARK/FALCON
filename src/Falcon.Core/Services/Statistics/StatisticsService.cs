using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.Statistics;

public sealed class StatisticsService : IStatisticsService
{
    public (double[][] NormalizedModels, double[] MassFactors) NormaliseSpec(ReadOnlySpan<double> dataFlux, double[][] modelFlux)
    {
        double dataNorm = Median(dataFlux);
        int numMods = modelFlux.Length;
        var normed = new double[numMods][];
        var massFactors = new double[numMods];
        for (int m = 0; m < numMods; m++)
        {
            double modelNorm = Median(modelFlux[m]);
            massFactors[m] = dataNorm / modelNorm;
            normed[m] = modelFlux[m].Select(v => v / modelNorm * dataNorm).ToArray();
        }
        return (normed, massFactors);
    }

    public (double[] MatchedWave, double[] MatchedData, double[] MatchedError, double[][] MatchedModels) MatchDataModels(
        double[] dataWave, double[] dataFlux, bool[] dataFlags, double[] errorFlux,
        double[] modelWave, double[][] modelFlux, double minWave, double maxWave)
    {
        int numModels = modelFlux.Length;
        double minW = Math.Max(dataWave.Where((_, i) => dataFlags[i]).DefaultIfEmpty(dataWave.Min()).Min(), Math.Max(modelWave.Min(), minWave));
        double maxW = Math.Min(dataWave.Where((_, i) => dataFlags[i]).DefaultIfEmpty(dataWave.Max()).Max(), Math.Min(modelWave.Max(), maxWave));

        var locModel = modelWave.Select(w => w >= minW && w <= maxW).ToArray();
        int numMod = locModel.Count(x => x);
        if (numMod == 0) throw new ArgumentException("Wavelength range outside model coverage.");
        var trimmedModelWave = modelWave.Where((w, i) => locModel[i]).ToArray();
        var trimmedModelFlux = new double[numModels][];
        for (int m = 0; m < numModels; m++)
            trimmedModelFlux[m] = modelFlux[m].Where((_, i) => locModel[i]).ToArray();

        var locData = dataWave.Select(w => w >= minW && w <= maxW).ToArray();
        int numDat = locData.Count(x => x);
        if (numDat == 0) throw new ArgumentException("Wavelength range outside data coverage.");
        var trimmedDataWave = dataWave.Where((_, i) => locData[i]).ToArray();
        var trimmedDataFlux = dataFlux.Where((_, i) => locData[i]).ToArray();
        var trimmedError = errorFlux.Where((_, i) => locData[i]).ToArray();

        if (numMod >= numDat)
        {
            var bisect = BisectArray(trimmedDataWave);
            var matchedModel = new double[numModels][];
            for (int m = 0; m < numModels; m++)
                matchedModel[m] = IntegrateToBins(trimmedModelWave, trimmedModelFlux[m], bisect);
            return (trimmedDataWave[1..^1], trimmedDataFlux[1..^1], trimmedError[1..^1], matchedModel);
        }

        var bisectModel = BisectArray(trimmedModelWave);
        var matchedData = IntegrateToBins(trimmedDataWave, trimmedDataFlux, bisectModel);
        var matchedError = IntegrateToBins(trimmedDataWave, trimmedError, bisectModel);
        var matchedModels = new double[numModels][];
        for (int m = 0; m < numModels; m++)
            matchedModels[m] = trimmedModelFlux[m][1..^1];
        return (trimmedModelWave[1..^1], matchedData, matchedError, matchedModels);
    }

    public double[] ConvertChisToProbs(ReadOnlySpan<double> chiSquared)
    {
        if (chiSquared.IsEmpty) return Array.Empty<double>();
        double minChi = chiSquared.ToArray().Min();
        var scaled = chiSquared.ToArray().Select(c => Math.Exp(-0.5 * (c - minChi))).ToArray();
        double sum = scaled.Sum();
        return scaled.Select(v => v / sum).ToArray();
    }

    public Dictionary<string, double> CalculateAveragesPdf(double[] probs, double[] lightWeights, double[] massWeights, double[] age, double[] metal, double fluxUnits, double distLumCm, double[] unnormMass, int sampling)
    {
        var result = new Dictionary<string, double>();
        result["light_age"] = Dot(probs, lightWeights, age);
        result["light_metal"] = Dot(probs, lightWeights, metal);
        result["mass_age"] = Dot(probs, massWeights, age);
        result["mass_metal"] = Dot(probs, massWeights, metal);
        double conversion = fluxUnits * 4 * Math.PI * distLumCm * distLumCm;
        double totMass = unnormMass.Zip(probs, (m, p) => m * p).Sum() * conversion;
        result["stellar_mass"] = totMass;
        return result;
    }

    public double[] LightWeightsToMass(double[] lightWeights, double[] massFactors)
    {
        return lightWeights.Zip(massFactors, (w, f) => w * f).ToArray();
    }

    private static double Dot(double[] probs, double[] weights, double[] values)
    {
        double s = 0;
        for (int i = 0; i < probs.Length; i++) s += probs[i] * weights[i] * values[i];
        return s;
    }

    private static double Median(ReadOnlySpan<double> data)
    {
        var sorted = data.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2 : sorted[n / 2];
    }

    private static double[] BisectArray(double[] wave)
    {
        var bins = new List<double> { wave[0] };
        for (int i = 1; i < wave.Length; i++)
            bins.Add((wave[i] + wave[i - 1]) / 2);
        bins.Add(wave[^1]);
        return bins.ToArray();
    }

    private static double[] IntegrateToBins(double[] wave, double[] flux, double[] bisect)
    {
        var combinedWave = wave.Concat(bisect).OrderBy(x => x).ToArray();
        var interp = combinedWave.Select(x =>
        {
            int idx = Array.BinarySearch(wave, x);
            if (idx >= 0) return flux[idx];
            idx = ~idx;
            if (idx <= 0) return flux[0];
            if (idx >= wave.Length) return flux[^1];
            double t = (x - wave[idx - 1]) / (wave[idx] - wave[idx - 1]);
            return flux[idx - 1] * (1 - t) + flux[idx] * t;
        }).ToArray();

        var output = new double[bisect.Length - 1];
        for (int l = 0; l < output.Length; l++)
        {
            int start = Array.BinarySearch(combinedWave, bisect[l]);
            if (start < 0) start = ~start;
            int end = Array.BinarySearch(combinedWave, bisect[l + 1]);
            if (end < 0) end = ~end;
            if (end >= combinedWave.Length) { output[l] = l > 0 ? output[l - 1] : 0; continue; }
            double integral = 0;
            for (int i = start; i < end; i++)
                integral += interp[i] * (combinedWave[i + 1] - combinedWave[i]);
            output[l] = integral / (bisect[l + 1] - bisect[l]);
        }
        return output;
    }
}
