using System.Text.Json;
using System.Text.Json.Serialization;
using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.IO;

public static class FalconResultsWriter
{
    public static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true,
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
        NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals
    };

    public static string ResultsJsonPath(string outputBasePath) =>
        Path.ChangeExtension(outputBasePath, ".results.json");

    public static void Write(FireflyFitOutput output, float redshift, string engine = "firefly")
    {
        Directory.CreateDirectory(Path.GetDirectoryName(output.ResultsPath)!);
        var doc = BuildDocument(output, redshift, engine);
        File.WriteAllText(output.ResultsPath, JsonSerializer.Serialize(doc, JsonOptions));
    }

    public static FalconResultsDocument BuildDocument(FireflyFitOutput output, float redshift, string engine)
    {
        int dof = output.MatchedWave.Length;
        var (waveMin, waveMax) = ComputeModelWaveRange(output.MatchedWave, output.BestFitSpectrum);
        return new FalconResultsDocument
        {
            SourceFile = output.SourceFile,
            Engine = engine,
            Redshift = redshift,
            Ebv = output.BestEbv,
            Converged = output.Converged,
            DegreesOfFreedom = dof,
            MinChiSquared = output.MinChiSquared,
            WavelengthAngstrom = output.MatchedWave,
            DataFlux = output.MatchedData,
            ErrorFlux = output.MatchedError,
            ModelWaveMin = waveMin,
            ModelWaveMax = waveMax,
            PdfAverages = output.Averages,
            TopSolutions = output.TopSolutions
        };
    }

    public static (double Min, double Max) ComputeModelWaveRange(double[] wave, double[] model)
    {
        if (wave.Length == 0 || model.Length == 0)
            return (0, 0);

        int n = Math.Min(wave.Length, model.Length);
        double maxAbs = 0;
        for (int i = 0; i < n; i++)
        {
            double a = Math.Abs(model[i]);
            if (a > maxAbs) maxAbs = a;
        }

        double eps = maxAbs * 1e-12;
        int iMin = 0;
        int iMax = n - 1;
        while (iMin < iMax && Math.Abs(model[iMin]) <= eps) iMin++;
        while (iMax > iMin && Math.Abs(model[iMax]) <= eps) iMax--;
        return (wave[iMin], wave[iMax]);
    }
}

public static class FalconResultsReader
{
    public static FalconResultsDocument Load(string jsonPath)
    {
        return JsonSerializer.Deserialize<FalconResultsDocument>(
            File.ReadAllText(jsonPath), FalconResultsWriter.JsonOptions)
            ?? throw new InvalidDataException($"Could not parse results: {jsonPath}");
    }

    public static IEnumerable<string> FindResultFiles(string outputDirectory)
    {
        if (!Directory.Exists(outputDirectory)) yield break;
        foreach (var f in Directory.EnumerateFiles(outputDirectory, "*.results.json", SearchOption.TopDirectoryOnly))
            yield return f;
    }
}
