namespace Falcon.UI.Services;

public enum PlotRangeMode
{
    /// <summary>Show full data wavelength range (model may be zero at edges).</summary>
    Full,
    /// <summary>Crop to wavelength range where model flux is non-zero (dissertation style).</summary>
    SnipToModelOverlap
}
