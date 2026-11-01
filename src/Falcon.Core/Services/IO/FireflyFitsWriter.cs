using System.Text.Json;
using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.IO;

/// <summary>Writes firefly-compatible header sidecar; full spectra live in .results.json.</summary>
public sealed class FireflyFitsWriter : IFireflyFitsWriter
{
    public void WriteResults(string outputPath, FireflyFitOutput output, float redshift, MassLossHeaderValues? massLoss = null)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        var headers = BuildHeaders(output, redshift, massLoss, converged: output.Converged);
        string jsonPath = Path.ChangeExtension(outputPath, ".json");
        File.WriteAllText(jsonPath, JsonSerializer.Serialize(headers, new JsonSerializerOptions { WriteIndented = true }));
        File.WriteAllText(outputPath + ".headers.txt", string.Join('\n', headers.Select(kv => $"HIERARCH {kv.Key} = {kv.Value}")));
    }

    public void WriteDummyResults(string outputPath, float redshift)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        const double defaultValue = 999.0;
        var headers = new Dictionary<string, object>
        {
            ["redshift"] = redshift,
            ["converged"] = "False",
            ["age_lightW"] = defaultValue,
            ["metallicity_lightW"] = defaultValue,
            ["total_mass"] = defaultValue,
            ["EBV"] = defaultValue
        };
        File.WriteAllText(Path.ChangeExtension(outputPath, ".json"), JsonSerializer.Serialize(headers));
    }

    private static Dictionary<string, object> BuildHeaders(
        FireflyFitOutput output, float redshift, MassLossHeaderValues? massLoss, bool converged)
    {
        var headers = new Dictionary<string, object>
        {
            ["redshift"] = redshift,
            ["converged"] = converged ? "True" : "False",
            ["EBV"] = output.BestEbv,
            ["min_chi_squared"] = output.MinChiSquared,
            ["top_solutions_count"] = output.TopSolutions.Count
        };
        foreach (var (key, value) in output.Averages)
            headers[key] = value;
        if (massLoss is not null)
        {
            headers["stellar_mass_ml"] = massLoss.CombinedMlAlive + massLoss.CombinedMlWd + massLoss.CombinedMlNs + massLoss.CombinedMlBh;
            headers["living_stars_mass"] = massLoss.CombinedMlAlive;
            headers["remnant_mass"] = massLoss.CombinedMlWd + massLoss.CombinedMlNs + massLoss.CombinedMlBh;
            headers["mass_of_ejecta"] = massLoss.CombinedGasFraction;
        }
        return headers;
    }
}
