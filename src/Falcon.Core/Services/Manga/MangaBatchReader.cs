using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.Manga;

/// <summary>MaNGA cube batch reader — post core-parity milestone.</summary>
public interface IMangaBatchReader
{
    IEnumerable<string> EnumerateSpectraFromCube(string cubePath);
}

public sealed class MangaBatchReader : IMangaBatchReader
{
    public IEnumerable<string> EnumerateSpectraFromCube(string cubePath)
    {
        if (!File.Exists(cubePath))
            throw new FileNotFoundException("MaNGA cube not found", cubePath);
        // Placeholder: extract per-spaxel spectra to temp FITS for pipeline consumption.
        yield break;
    }
}
