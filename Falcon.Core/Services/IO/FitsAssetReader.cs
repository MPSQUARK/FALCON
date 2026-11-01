using CFITS;
using System.IO.Compression;

namespace Falcon.Core.Services.IO;

internal static class FitsAssetReader
{
    public static double[] ReadImage(string path, int hduIndex = 0)
    {
        if (hduIndex == 0)
        {
            DustMapFile file = FitsGeneratedReader.Read<DustMapFile>(path);
            float[,] pixels = file.Image.Pixels
                ?? throw new InvalidDataException("Dust-map image has no pixel data.");
            var values = new double[pixels.Length];
            for (int row = 0; row < pixels.GetLength(0); row++)
                for (int column = 0; column < pixels.GetLength(1); column++)
                    values[row * pixels.GetLength(1) + column] = pixels[row, column];
            return values;
        }
        return ReadImages(path, hduIndex)[hduIndex];
    }

    public static IReadOnlyDictionary<int, double[]> ReadImages(string path, params int[] hduIndices)
    {
        using FitsDocument document = OpenDocument(path);
        var images = new Dictionary<int, double[]>(hduIndices.Length);
        foreach (int hduIndex in hduIndices)
        {
            FitsImage image = document.Hdus[hduIndex].Image
                ?? throw new InvalidDataException($"FITS HDU {hduIndex} does not contain an image.");
            images[hduIndex] = image.ReadPhysical();
        }
        return images;
    }

    private static FitsDocument OpenDocument(string path)
    {
        if (!path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
            return FitsDocument.Open(path);

        return FitsDocument.Parse(FitsInput.ReadGzip(path));
    }
}

internal static class FitsGeneratedReader
{
    public static T Read<T>(string path) where T : CFITS.IFitsDocument<T>
    {
        if (!path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
            using (FileStream stream = File.OpenRead(path))
                return CFITS.FitsParser.Deserialize<T>(stream);

        using var decompressed = new MemoryStream(FitsInput.ReadGzip(path), writable: false);
        return CFITS.FitsParser.Deserialize<T>(decompressed);
    }
}

internal static class FitsInput
{
    public static byte[] ReadGzip(string path)
    {
        using FileStream input = File.OpenRead(path);
        using var gzip = new GZipStream(input, CompressionMode.Decompress);
        using var output = new MemoryStream();
        gzip.CopyTo(output);
        return output.ToArray();
    }
}
