using System.Text.Json;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Golden;

/// <summary>
/// Sanity checks against FIREFLY results from the thesis reference table.
/// Full e2e parity runs are optional (require MaStar HDF5 + GPU).
/// </summary>
public class ThesisFireflyReferenceTests
{
    private static JsonDocument LoadThesisRef()
    {
        string path = Path.Combine(AppContext.BaseDirectory, "TestData", "golden_thesis_firefly.json");
        return JsonDocument.Parse(File.ReadAllText(path));
    }

    [Fact]
    public void ThesisReference_ContainsSixGalaxies()
    {
        using var doc = LoadThesisRef();
        doc.RootElement.GetProperty("galaxies").GetArrayLength().Should().Be(6);
    }

    [Theory]
    [InlineData("0266-51630-0016", 1.06, 0.349)]
    [InlineData("0266-51630-0022", 9.96, 0.233)]
    [InlineData("0266-51630-0169", 11.98, 0.075)]
    public void ThesisReference_FireflyAgesAreInGyrRange(string id, double ageGyr, double metalDex)
    {
        using var doc = LoadThesisRef();
        var galaxy = doc.RootElement.GetProperty("galaxies").EnumerateArray()
            .First(g => g.GetProperty("id").GetString() == id);
        galaxy.GetProperty("firefly_age_gyr").GetDouble().Should().BeApproximately(ageGyr, 0.01);
        galaxy.GetProperty("firefly_metal_dex").GetDouble().Should().BeApproximately(metalDex, 0.01);
        ageGyr.Should().BeGreaterThan(0.1);
        ageGyr.Should().BeLessThan(14);
    }

    [Fact]
    public void ResultsJson_ModelFlux_IsNonZeroWhenPresent()
    {
        string? resultsPath = FindResultsJson("spec-0266-51630-0016.results.json");
        if (resultsPath is null) return;

        var doc = FalconResultsReader.Load(resultsPath);
        var top = doc.TopSolutions.FirstOrDefault();
        top.Should().NotBeNull();
        top!.ModelFlux.Should().NotBeEmpty();
        top.ModelFlux.Max(Math.Abs).Should().BeGreaterThan(0,
            "model flux must be non-zero for plotting; re-run fit after BuildFireflyContext fix");
    }

    [Fact]
    public void ResultsJson_LightAge_IsReasonableGyr_WhenPresent()
    {
        string? resultsPath = FindResultsJson("spec-0266-51630-0016.results.json");
        if (resultsPath is null) return;

        var doc = FalconResultsReader.Load(resultsPath);
        doc.PdfAverages.TryGetValue("light_age", out double lightAge).Should().BeTrue();
        lightAge.Should().BeGreaterThan(0.01, "ages are in Gyr, not years");
        lightAge.Should().BeLessThan(14, "ages should not exceed cosmic age in Gyr");
    }

    private static string? FindResultsJson(string fileName)
    {
        string[] candidates =
        [
            Path.Combine(AppContext.BaseDirectory, "output", fileName),
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "src", "Falcon.Cli", "bin", "Debug", "net10.0", "output", fileName)),
        ];
        return candidates.FirstOrDefault(File.Exists);
    }
}
