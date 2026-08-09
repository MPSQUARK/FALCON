namespace Falcon.Core.Configuration;

public sealed class FireflyOptions
{
    public InputOptions Input { get; set; } = new();
    public ModelOptions Model { get; set; } = new();
    public FittingOptions Fitting { get; set; } = new();
    public DustOptions Dust { get; set; } = new();
    public MaskingOptions Masking { get; set; } = new();
    public PathOptions Paths { get; set; } = new();
}

public sealed class InputOptions
{
    public string SpectrumPath { get; set; } = "assets/sample-data";
    public string InputGlob { get; set; } = "*.fits";
    public float Redshift { get; set; } = 0.021275453f;
    public float Ra { get; set; } = 145.89219f;
    public float Dec { get; set; } = 0.059372f;
    public float VelocityDispersion { get; set; } = 135.89957f;
    public float InstrumentResolution { get; set; } = 2000f;
    public sbyte FluxScaleFactor { get; set; } = -17;
    public bool DataMediumAir { get; set; }
    public bool WriteResults { get; set; } = true;
}

public sealed class ModelOptions
{
    public string Library { get; set; } = "MaStar-E";
    public string Imf { get; set; } = "Kroupa";
    public float AgeMin { get; set; }
    public float AgeMax { get; set; } = 14.7e9f;
    public float MetallicityMin { get; set; } = -3f;
    public float MetallicityMax { get; set; } = 3f;
    public bool DowngradeModels { get; set; } = true;
    public byte ModelKey { get; set; } = 0b00010_010;
}

public sealed class FittingOptions
{
    public int MaxIterations { get; set; } = 10;
    public int PdfSampling { get; set; } = 200;
    public float MaxEbv { get; set; } = 1.0f;
    public int NumDustVals { get; set; } = 200;
    public float DustSmoothingLength { get; set; } = 100f;
}

public sealed class DustOptions
{
    public bool MilkyWayReddening { get; set; } = true;
    public string HpfMode { get; set; } = "on";
    public string DustLaw { get; set; } = "calzetti";
}

public sealed class MaskingOptions
{
    public ushort MaskedAngstroms { get; set; } = 20;
    public float[] EmissionLines { get; set; } =
    [
        6562.80f, 6583.34f, 6547.96f, 4861.32f, 5006.77f, 4958.83f
    ];
}

public sealed class PathOptions
{
    public string StellarPopulationModels { get; set; } = "assets/stellar-population-models";
    public string DustMaps { get; set; } = "assets/dust-maps";
    public string Output { get; set; } = "output";
    public string MassLossTables { get; set; } = "assets/mass-loss";
}
