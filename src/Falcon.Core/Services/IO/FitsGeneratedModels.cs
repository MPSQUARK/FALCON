using CFITS;

namespace Falcon.Core.Services.IO;

[FitsDocument]
internal sealed partial class SdssSpectrumFile
{
    [FitsPrimaryHeader]
    public SdssPrimaryHeader Primary { get; set; } = null!;

    [FitsTable(1)]
    [FitsRequired]
    public SdssSpectrumTable Spectrum { get; set; } = null!;

    [FitsTable(2)]
    public SdssSpecObjTable SpecObj { get; set; } = new();
}

internal sealed class SdssPrimaryHeader
{
    [FitsKeyword("RA")]
    public double Ra { get; set; }

    [FitsKeyword("DEC")]
    public double Dec { get; set; }
}

internal sealed class SdssSpectrumTable
{
    public IReadOnlyList<SdssSpectrumRow> Rows { get; set; } = [];
}

internal sealed class SdssSpectrumRow
{
    [FitsColumn("flux")]
    [FitsRequired]
    public float Flux { get; set; }

    [FitsColumn("loglam")]
    [FitsRequired]
    public float LogWavelength { get; set; }

    [FitsColumn("ivar")]
    [FitsRequired]
    public float InverseVariance { get; set; }
}

internal sealed class SdssSpecObjTable
{
    public IReadOnlyList<SdssSpecObjRow> Rows { get; set; } = [];
}

internal sealed class SdssSpecObjRow
{
    [FitsColumn("Z")]
    public float Redshift { get; set; }

    [FitsColumn("Z_NOQSO")]
    public float RedshiftNoQso { get; set; }

    [FitsColumn("VDISP")]
    public float VelocityDispersion { get; set; }
}

[FitsDocument]
internal sealed partial class DustMapFile
{
    [FitsImage]
    [FitsRequired]
    public DustMapImage Image { get; set; } = null!;
}

internal sealed class DustMapImage
{
    [FitsData]
    [FitsRequired]
    public float[,]? Pixels { get; set; }
}
