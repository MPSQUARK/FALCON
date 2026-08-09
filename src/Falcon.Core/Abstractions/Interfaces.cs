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
    public int FitPerIterationCap { get; init; }
    public int MixTopCount { get; init; } = 10;
    public double SigmaClipLow { get; init; } = 3.0;
    public double SigmaClipHigh { get; init; } = 3.0;
    public IReadOnlyList<double> MixQuantities { get; init; } =
        [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0];
}

public interface IChiSquaredEvaluator
{
    ReadOnlyMemory<double> ChiModels { get; }
    int NumModels { get; }
    int NumWave { get; }
    double Evaluate(ReadOnlySpan<double> weights, int branchNumber, Span<bool> clippedMask);
}

public interface IDustService
{
    double[] Hpf(ReadOnlySpan<double> flux, int windowSize = 20, int wStart = 20);
    double[] Unred(ReadOnlySpan<double> wavelengthAngstrom, double ebv, string law = "calzetti");
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
    double[] ConvertChisToProbs(ReadOnlySpan<double> chiSquared);
    Dictionary<string, double> CalculateAveragesPdf(double[] probs, double[] lightWeights, double[] massWeights, double[] age, double[] metal, double fluxUnits, double distLumCm, double[] unnormMass, int sampling);
    double[] LightWeightsToMass(double[] lightWeights, double[] massFactors);
}

public interface IResolutionService
{
    // Marker for instrument downgrade — implemented by ResolutionService in Services.Instrument
}

public interface ISpectrumReader
{
    (Vector Wavelength, Vector Flux, Vector Error, float Redshift, float Vdisp, float Ra, float Dec) ReadFits(GPU gpu, string path);
}

public interface ISpectralFitPipeline
{
    void FitFile(string fitsPath);
    void FitBatch(IEnumerable<string> fitsPaths, int maxParallelism = -1);
}
