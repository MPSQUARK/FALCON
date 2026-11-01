using System.Text.Json;
using Falcon.Core.Abstractions;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class FireflyFitsWriterTests
{
    [Fact]
    public void WriteResults_CreatesJsonSidecar()
    {
        var dir = CreateTempDir();
        var path = Path.Combine(dir, "result.fits");
        var output = SampleOutput();
        new FireflyFitsWriter().WriteResults(path, output, redshift: 0.02f);
        File.Exists(Path.ChangeExtension(path, ".json")).Should().BeTrue();
    }

    [Fact]
    public void WriteResults_JsonContainsExpectedHeaders()
    {
        var dir = CreateTempDir();
        var path = Path.Combine(dir, "result.fits");
        new FireflyFitsWriter().WriteResults(path, SampleOutput(), 0.02f);
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.ChangeExtension(path, ".json")));
        doc.RootElement.GetProperty("redshift").GetDouble().Should().BeApproximately(0.02, 1e-6);
        doc.RootElement.GetProperty("EBV").GetDouble().Should().BeApproximately(0.1, 1e-6);
        doc.RootElement.GetProperty("converged").GetString().Should().Be("True");
    }

    [Fact]
    public void WriteResults_WithMassLoss_AddsMassHeaders()
    {
        var dir = CreateTempDir();
        var path = Path.Combine(dir, "result.fits");
        var massLoss = new MassLossHeaderValues
        {
            CombinedMlAlive = 1e10,
            CombinedMlWd = 1e9,
            CombinedMlNs = 1e8,
            CombinedMlBh = 1e7,
            CombinedGasFraction = 1e6
        };
        new FireflyFitsWriter().WriteResults(path, SampleOutput(), 0.02f, massLoss);
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.ChangeExtension(path, ".json")));
        doc.RootElement.GetProperty("living_stars_mass").GetDouble().Should().Be(1e10);
        doc.RootElement.GetProperty("stellar_mass_ml").GetDouble().Should().BeApproximately(1e10 + 1e9 + 1e8 + 1e7, 1e5);
    }

    [Fact]
    public void WriteDummyResults_WritesSentinelValues()
    {
        var dir = CreateTempDir();
        var path = Path.Combine(dir, "dummy.fits");
        new FireflyFitsWriter().WriteDummyResults(path, 0.05f);
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.ChangeExtension(path, ".json")));
        doc.RootElement.GetProperty("age_lightW").GetDouble().Should().Be(999.0);
        doc.RootElement.GetProperty("converged").GetString().Should().Be("False");
    }

    [Fact]
    public void WriteResults_CreatesHeadersTextFile()
    {
        var dir = CreateTempDir();
        var path = Path.Combine(dir, "result.fits");
        new FireflyFitsWriter().WriteResults(path, SampleOutput(), 0.02f);
        File.Exists(path + ".headers.txt").Should().BeTrue();
        File.ReadAllText(path + ".headers.txt").Should().Contain("HIERARCH");
    }

    private static FireflyFitOutput SampleOutput() => new()
    {
        SourceFile = "spec.fits",
        Converged = true,
        BestEbv = 0.1,
        MinChiSquared = 1.23,
        BestLightWeights = [0.5, 0.5],
        BestMassWeights = [0.5, 0.5],
        Averages = new Dictionary<string, double> { ["light_age"] = 5e9, ["stellar_mass"] = 1e10 },
        MatchedWave = [4000, 4100],
        MatchedData = [1, 2],
        MatchedError = [0.1, 0.1],
        BestFitSpectrum = [1, 2],
        OutputPath = "out",
        ResultsPath = "out.results.json",
        TopSolutions = []
    };

    private static string CreateTempDir()
    {
        var dir = Path.Combine(Path.GetTempPath(), "falcon-writer-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        return dir;
    }
}
