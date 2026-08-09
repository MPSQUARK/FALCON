using FluentAssertions;

namespace Falcon.Tests.Golden;

public class GoldenManifestTests
{
    [Fact]
    public void GoldenManifest_Exists()
    {
        var path = Path.Combine(AppContext.BaseDirectory, "TestData", "golden_manifest.json");
        File.Exists(path).Should().BeTrue("bundled golden manifest should ship with tests");
    }
}
