using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.IO;

/// <summary>Canonical JSON results format for plotting and UI display.</summary>
public sealed class FalconResultsDocument
{
    public int Version { get; set; } = 1;
    public required string SourceFile { get; set; }
    public required string Engine { get; set; }
    public float Redshift { get; set; }
    public double Ebv { get; set; }
    public bool Converged { get; set; }
    public int DegreesOfFreedom { get; set; }
    public double MinChiSquared { get; set; }
    public double[] WavelengthAngstrom { get; set; } = [];
    public double[] DataFlux { get; set; } = [];
    public double[] ErrorFlux { get; set; } = [];
    /// <summary>Rest wavelength where model flux becomes non-negligible (for plot snip).</summary>
    public double ModelWaveMin { get; set; }
    /// <summary>Rest wavelength where model flux becomes non-negligible (for plot snip).</summary>
    public double ModelWaveMax { get; set; }
    public Dictionary<string, double> PdfAverages { get; set; } = new();
    public List<FalconSolutionRecord> TopSolutions { get; set; } = [];
}
