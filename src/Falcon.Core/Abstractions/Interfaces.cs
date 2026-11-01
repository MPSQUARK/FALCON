namespace Falcon.Core.Abstractions;

public interface ISpectralFitter
{
    FitResult Fit(SpectralFitInput input, FitterOptions? options = null);
}

public sealed class SpectralFitInput
{
    public required ReadOnlyMemory<double> Wavelength { get; init; }
    public required ReadOnlyMemory<double> Data { get; init; }
    public required ReadOnlyMemory<double> Error { get; init; }
    public required ReadOnlyMemory<double> Models { get; init; }

    public int NumModels => Data.Length > 0 && Models.Length % Data.Length == 0
        ? Models.Length / Data.Length
        : throw new InvalidOperationException("Models length must be an integer multiple of data length.");

    public int NumWave => Data.Length;

    /// <summary>Build fit input from float32 spectra (promoted to double only at property boundaries).</summary>
    public static SpectralFitInput FromFloat(
        ReadOnlyMemory<float> wavelength,
        ReadOnlyMemory<float> data,
        ReadOnlyMemory<float> error,
        ReadOnlyMemory<float> models)
    {
        return new SpectralFitInput
        {
            Wavelength = ToDoubleMemory(wavelength),
            Data = ToDoubleMemory(data),
            Error = ToDoubleMemory(error),
            Models = ToDoubleMemory(models),
        };
    }

    static double[] ToDoubleMemory(ReadOnlyMemory<float> values)
    {
        var src = values.Span;
        var dst = new double[src.Length];
        for (int i = 0; i < src.Length; i++)
            dst[i] = src[i];
        return dst;
    }

    public void Validate()
    {
        if (Wavelength.Length != NumWave)
            throw new ArgumentException("Wavelength length must match data length.");
        if (Error.Length != NumWave)
            throw new ArgumentException("Error length must match data length.");
        if (Models.Length != NumModels * NumWave)
            throw new ArgumentException("Models length must equal NumModels * NumWave.");
    }
}

public sealed class FitResult
{
    public required IReadOnlyList<FitSolution> Solutions { get; init; }
    public required double MinChiSquared { get; init; }
    public required bool[] ClippedMask { get; init; }
    public required int IterationCount { get; init; }
    public required bool Converged { get; init; }
}

public sealed class FitSolution
{
    public required int Index { get; init; }
    public required int BranchNumber { get; init; }
    public required double[] Weights { get; init; }
    public required double ChiSquared { get; init; }
    public bool[]? ClippedMask { get; init; }
}

public sealed class FitterOptions
{
    public int MaxIterations { get; init; } = 10;
    public int FitPerIterationCap { get; init; } = 1000;
    public int MixTopCount { get; init; } = 10;
    /// <summary>Maximum distinct SSP components per blend (firefly explores up to five).</summary>
    public int MaxSspComponents { get; init; } = 5;
    public double SigmaClipLow { get; init; } = 3.0;
    public double SigmaClipHigh { get; init; } = 3.0;
    public IReadOnlyList<double> MixQuantities { get; init; } =
        [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0];
    /// <summary>Use GPU batched χ² blend when a GPU is available (default false — firefly-equivalent CPU dot is faster for typical Mix eval counts).</summary>
    public bool UseGpuEvaluator { get; init; } = false;
    /// <summary>Replace firefly mix q-loop with golden-section search for two-active-component blends only.</summary>
    public bool ExperimentalMixSearch { get; init; } = false;
    public int ExperimentalMixSearchMaxEvals { get; init; } = 64;
}

public interface IChiSquaredEvaluator : IDisposable
{
    ReadOnlyMemory<double> ChiModels { get; }
    int NumModels { get; }
    int NumWave { get; }
    bool SupportsBatchEvaluation { get; }
    bool[] SharedClippedMask { get; }
    void SetSharedClippedMask(ReadOnlySpan<bool> mask);
    double Evaluate(ReadOnlySpan<double> weights, int branchNumber, Span<bool> clippedMask);

    /// <summary>Score multiple weight vectors; returns number of entries written to <paramref name="chiSquaredOut"/>.</summary>
    int EvaluateBatch(
        ReadOnlySpan<double> flatWeights,
        int batchSize,
        ReadOnlySpan<int> branchNumbers,
        Span<double> chiSquaredOut);
}

public interface IDustService
{
    double[] Hpf(ReadOnlySpan<double> flux, int windowSize = 20, int wStart = 20);
    double[] Unred(ReadOnlySpan<double> wavelengthAngstrom, double ebv, string law = "calzetti");
    /// <summary>Fitzpatrick (1999) MW dereddening — firefly <c>unred()</c> for model correction.</summary>
    double[] FitzpatrickUnred(ReadOnlySpan<double> wavelengthAngstrom, double ebv, double rV = 3.1);
    double GetSfdEbv(double longitudeDeg, double latitudeDeg, float[,] ngp, float[,] sgp);
    DustAttenuationResult DetermineAttenuation(DustAttenuationInput input, ISpectralFitter fitter);
}

public sealed class DustAttenuationInput
{
    public required double[] Wavelength { get; init; }
    public required double[] DataFlux { get; init; }
    public required double[] ErrorFlux { get; init; }
    public required double[][] ModelFlux { get; init; }
    public required double[] Age { get; init; }
    public required double[] Metal { get; init; }
    public float MaxEbv { get; init; } = 1f;
    public int NumDustVals { get; init; } = 200;
    public float SmoothingLength { get; init; } = 100f;
    public string DustLaw { get; init; } = "calzetti";
}

public sealed class DustAttenuationResult
{
    public required double Ebv { get; init; }
    public required double[] AttenuationCurve { get; init; }
    public required double AvAgeHpf { get; init; }
    public required double AvMetalHpf { get; init; }
}

public interface IStatisticsService
{
    (double[][] NormalizedModels, double[] MassFactors) NormaliseSpec(ReadOnlySpan<double> dataFlux, double[][] modelFlux);
    (double[] MatchedWave, double[] MatchedData, double[] MatchedError, double[][] MatchedModels) MatchDataModels(
        double[] dataWave, double[] dataFlux, bool[] dataFlags, double[] errorFlux,
        double[] modelWave, double[][] modelFlux, double minWave, double maxWave);
    double[] ConvertChisToProbs(ReadOnlySpan<double> chiSquared, int dof);
    Dictionary<string, double> CalculateAveragesPdf(
        double[] probs, double[][] lightWeights, double[][] massWeights, double[][] unnormMass,
        double[] age, double[] metal, double fluxUnits, double distLumCm, int sampling);
    (double[][] UnnormMass, double[][] MassWeights) LightWeightsToMass(double[][] lightWeights, double[] massFactors);
    double[] LightWeightsToMassRow(double[] lightWeights, double[] massFactors);
}

public sealed class FireflyFitOutput
{
    public required string SourceFile { get; init; }
    public required bool Converged { get; init; }
    public required double BestEbv { get; init; }
    public required double MinChiSquared { get; init; }
    public required double[] BestLightWeights { get; init; }
    public required double[] BestMassWeights { get; init; }
    public required Dictionary<string, double> Averages { get; init; }
    public required double[] MatchedWave { get; init; }
    public required double[] MatchedData { get; init; }
    public required double[] MatchedError { get; init; }
    public required double[] BestFitSpectrum { get; init; }
    /// <summary>Base path for FITS/header sidecar (legacy).</summary>
    public required string OutputPath { get; init; }
    /// <summary>Path to <c>.results.json</c> with top-N solutions and spectra.</summary>
    public required string ResultsPath { get; init; }
    public required List<FalconSolutionRecord> TopSolutions { get; init; }
}

// Forward declare for Abstractions — implementation in Falcon.Core.Services.IO
public sealed class FalconSolutionRecord
{
    public int Rank { get; set; }
    public double ChiSquared { get; set; }
    public double ChiSquaredPerDof { get; set; }
    public double Probability { get; set; }
    public double LightWeightedAgeGyr { get; set; }
    public double LightWeightedMetal { get; set; }
    public double[] ModelFlux { get; set; } = [];
    public List<FalconSspComponent> Components { get; set; } = [];
}

public sealed class FalconSspComponent
{
    public int ModelIndex { get; set; }
    public double LightWeight { get; set; }
    public double AgeGyr { get; set; }
    public double Metal { get; set; }
}

public interface IFireflyFitsWriter
{
    void WriteResults(string outputPath, FireflyFitOutput output, float redshift, MassLossHeaderValues? massLoss = null);
    void WriteDummyResults(string outputPath, float redshift);
}

public sealed class MassLossHeaderValues
{
    public double CombinedMlAlive { get; init; }
    public double CombinedMlWd { get; init; }
    public double CombinedMlNs { get; init; }
    public double CombinedMlBh { get; init; }
    public double CombinedGasFraction { get; init; }
}

public interface IResolutionService
{
}

public interface ISpectrumReader
{
    (Vector Wavelength, Vector Flux, Vector Error, float Redshift, float Vdisp, float Ra, float Dec) ReadFits(GPU gpu, string path);
}

public interface ISpectralFitPipeline
{
    void FitFile(string fitsPath, CancellationToken cancellationToken = default);
    void FitBatch(IEnumerable<string> fitsPaths, CancellationToken cancellationToken = default, int maxParallelism = -1);
}
