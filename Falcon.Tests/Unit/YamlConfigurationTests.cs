using Falcon.Core.Configuration;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class YamlConfigurationTests
{
    private static string ExampleYamlPath()
    {
        var dir = AppContext.BaseDirectory;
        for (int i = 0; i < 8; i++)
        {
            var candidate = Path.GetFullPath(Path.Combine(dir, "..", "..", "..", "..", "..", "configs", "example.yaml"));
            if (File.Exists(candidate)) return candidate;
            dir = Path.GetDirectoryName(dir)!;
        }
        throw new FileNotFoundException("configs/example.yaml not found relative to test output");
    }

    [Fact]
    public void LoadFireflyYaml_LoadsExampleConfig()
    {
        var options = YamlConfigurationExtensions.LoadFireflyYaml(ExampleYamlPath());
        options.Should().NotBeNull();
        options.Input.Redshift.Should().BeApproximately(0.021275453f, 1e-6f);
    }

    [Fact]
    public void LoadFireflyYaml_MapsInputSection()
    {
        var o = YamlConfigurationExtensions.LoadFireflyYaml(ExampleYamlPath());
        o.Input.Ra.Should().BeApproximately(145.89219f, 1e-4f);
        o.Input.Dec.Should().BeApproximately(0.059372f, 1e-5f);
        o.Input.FluxScaleFactor.Should().Be(-17);
        o.Input.WriteResults.Should().BeTrue();
    }

    [Fact]
    public void LoadFireflyYaml_MapsModelAndFittingSections()
    {
        var o = YamlConfigurationExtensions.LoadFireflyYaml(ExampleYamlPath());
        o.Model.Library.Should().Be("MaStar-E");
        o.Model.Imf.Should().Be("Kroupa");
        o.Model.ModelKey.Should().Be(18);
        o.Fitting.MaxIterations.Should().Be(10);
        o.Fitting.MaxSspComponents.Should().Be(5);
        o.Fitting.PdfSampling.Should().Be(200);
    }

    [Fact]
    public void LoadFireflyYaml_MapsDustAndPaths()
    {
        var o = YamlConfigurationExtensions.LoadFireflyYaml(ExampleYamlPath());
        o.Dust.MilkyWayReddening.Should().BeTrue();
        o.Dust.HpfMode.Should().Be("on");
        o.Dust.DustLaw.Should().Be("calzetti");
        o.Paths.Output.Should().Be("output");
    }

    [Fact]
    public void LoadFireflyYaml_ExampleValuesMatchLegacyMappingExpectations()
    {
        var o = YamlConfigurationExtensions.LoadFireflyYaml(ExampleYamlPath());
        o.Dust.HpfMode.ToLowerInvariant().Should().Be("on");
        o.Model.Imf.Should().Be("Kroupa");
        o.Fitting.WaveLimits.Should().HaveCount(2);
        o.Fitting.WaveLimits[1].Should().BeApproximately(99_999_990, 1);
    }

    [Theory]
    [InlineData("allen", 1)]
    [InlineData("prevot", 2)]
    [InlineData("calzetti", 0)]
    public void YamlDustLaw_MapsToLegacyEncoding(string law, byte expected)
    {
        var options = new FireflyOptions { Dust = new DustOptions { DustLaw = law } };
        MapDustLaw(options.Dust.DustLaw).Should().Be(expected);
    }

    private static byte MapDustLaw(string law) => law.ToLowerInvariant() switch
    {
        "allen" => 1,
        "prevot" => 2,
        _ => 0
    };
}
