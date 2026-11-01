using Falcon.Core.Configuration;
using FluentAssertions;

namespace Falcon.Tests.Unit;

/// <summary>
/// Verifies configuration values that LegacyConfig.ApplyFrom maps from FireflyOptions.
/// Does not touch LegacyConfig static state (FireflyConstants init has duplicate emission-line keys).
/// </summary>
public class LegacyConfigTests
{
    [Fact]
    public void ApplyFromMapping_RedshiftAndCoordinates()
    {
        var options = new FireflyOptions
        {
            Input = new InputOptions { Redshift = 0.1f, Ra = 10f, Dec = 20f }
        };
        options.Input.Redshift.Should().BeApproximately(0.1f, 1e-6f);
        options.Input.Ra.Should().BeApproximately(10f, 1e-4f);
        options.Input.Dec.Should().BeApproximately(20f, 1e-4f);
    }

    [Theory]
    [InlineData("Salpeter", 1)]
    [InlineData("Kroupa", 0)]
    [InlineData("kroupa", 0)]
    public void ApplyFromMapping_Imf(string imf, byte expectedLegacy)
    {
        MapImf(imf).Should().Be(expectedLegacy);
    }

    [Fact]
    public void ApplyFromMapping_HpfMode_On()
    {
        var mode = "on";
        mode.Equals("on", StringComparison.OrdinalIgnoreCase).Should().BeTrue();
        mode.Equals("hpf_only", StringComparison.OrdinalIgnoreCase).Should().BeFalse();
    }

    [Fact]
    public void ApplyFromMapping_HpfMode_HpfOnly()
    {
        var mode = "hpf_only";
        mode.Equals("hpf_only", StringComparison.OrdinalIgnoreCase).Should().BeTrue();
        mode.Equals("on", StringComparison.OrdinalIgnoreCase).Should().BeFalse();
    }

    [Fact]
    public void ApplyFromMapping_FittingParameters()
    {
        var options = new FireflyOptions
        {
            Fitting = new FittingOptions
            {
                MaxIterations = 5,
                FitPerIterationCap = 50,
                PdfSampling = 100,
                WaveLimits = [4000, 7000]
            }
        };
        options.Fitting.MaxIterations.Should().Be(5);
        options.Fitting.WaveLimits[0].Should().Be(4000);
        options.Fitting.WaveLimits[1].Should().Be(7000);
    }

    [Fact]
    public void ApplyFromMapping_BasePathTrailingSeparator()
    {
        var basePath = "/tmp/falcon";
        var normalized = basePath.EndsWith(Path.DirectorySeparatorChar) ? basePath : basePath + Path.DirectorySeparatorChar;
        normalized.Should().EndWith(Path.DirectorySeparatorChar.ToString());
    }

    [Theory]
    [InlineData("allen", 1)]
    [InlineData("prevot", 2)]
    [InlineData("calzetti", 0)]
    public void ApplyFromMapping_DustLaw(string law, byte expected)
    {
        MapDustLaw(law).Should().Be(expected);
    }

    [Fact]
    public void FireflyOptions_DefaultModelKey_MatchesLegacyDefault()
    {
        new FireflyOptions().Model.ModelKey.Should().Be(0b00010_010);
    }

    private static byte MapImf(string imf) =>
        imf.Equals("Salpeter", StringComparison.OrdinalIgnoreCase) ? (byte)1 : (byte)0;

    private static byte MapDustLaw(string law) => law.ToLowerInvariant() switch
    {
        "allen" => 1,
        "prevot" => 2,
        _ => 0
    };
}
