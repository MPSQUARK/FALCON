using System.Text.Json;
using FluentAssertions;

namespace Falcon.Tests.Golden;

public class GoldenManifestTests
{
    [Fact]
    public void GoldenManifest_Exists()
    {
        ManifestPath().Should().NotBeNull();
        File.Exists(ManifestPath()).Should().BeTrue("bundled golden manifest should ship with tests");
    }

    [Fact]
    public void GoldenManifest_ListsAllGoldenFiles()
    {
        using var doc = JsonDocument.Parse(File.ReadAllText(ManifestPath()));
        var files = doc.RootElement.GetProperty("files").EnumerateArray().Select(e => e.GetString()!).ToArray();
        files.Should().NotBeEmpty();
        files.Should().Contain("golden_bisect.json");
        files.Should().Contain("golden_fitter.json");
    }

    [Fact]
    public void GoldenManifest_AllReferencedFilesExist()
    {
        using var doc = JsonDocument.Parse(File.ReadAllText(ManifestPath()));
        var testDataDir = Path.Combine(AppContext.BaseDirectory, "TestData");
        foreach (var file in doc.RootElement.GetProperty("files").EnumerateArray())
        {
            var name = file.GetString()!;
            File.Exists(Path.Combine(testDataDir, name)).Should().BeTrue($"golden file {name} should exist");
        }
    }

    [Fact]
    public void GoldenManifest_HasToleranceMetadata()
    {
        using var doc = JsonDocument.Parse(File.ReadAllText(ManifestPath()));
        doc.RootElement.GetProperty("tolerance_rtol").GetDouble().Should().BeGreaterThan(0);
        doc.RootElement.GetProperty("description").GetString().Should().NotBeNullOrWhiteSpace();
    }

    private static string ManifestPath() =>
        Path.Combine(AppContext.BaseDirectory, "TestData", "golden_manifest.json");
}
