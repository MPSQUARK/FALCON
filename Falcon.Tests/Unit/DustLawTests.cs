using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class DustLawTests : IDisposable
{
    private readonly GPU _gpu = FalconGpuHost.Create();
    private readonly DustService _dust;

    public DustLawTests() => _dust = new DustService(_gpu);

    public void Dispose() => _gpu.Dispose();
    private static readonly double[] Wave = [3500, 4500, 5500, 6500, 8500];

    [Theory]
    [InlineData("calzetti")]
    [InlineData("allen")]
    [InlineData("prevot")]
    public void Unred_AllLaws_ReturnPositiveAttenuation(string law)
    {
        var atten = _dust.Unred(Wave, ebv: 0.3, law);
        atten.Should().OnlyContain(v => v > 0 && v <= 1.0);
    }

    [Fact]
    public void Unred_Calzetti_HigherEbvMoreAttenuation()
    {
        var low = _dust.Unred(Wave, 0.1, "calzetti");
        var high = _dust.Unred(Wave, 0.5, "calzetti");
        for (int i = 0; i < Wave.Length; i++)
            high[i].Should().BeLessThan(low[i]);
    }

    [Fact]
    public void Unred_AllenDiffersFromCalzetti()
    {
        var cal = _dust.Unred(Wave, 0.2, "calzetti");
        var allen = _dust.Unred(Wave, 0.2, "allen");
        allen.Should().NotEqual(cal);
    }

    [Fact]
    public void Unred_PrevotDiffersFromCalzetti()
    {
        var cal = _dust.Unred(Wave, 0.2, "calzetti");
        var prev = _dust.Unred(Wave, 0.2, "prevot");
        prev.Should().NotEqual(cal);
    }

    [Theory]
    [InlineData("CALZETTI")]
    [InlineData("Allen")]
    [InlineData("PREVOT")]
    public void Unred_LawNameIsCaseInsensitive(string law)
    {
        var atten = _dust.Unred(Wave, 0.1, law);
        atten.Length.Should().Be(Wave.Length);
    }

    [Fact]
    public void Unred_ZeroEbv_ReturnsUnity()
    {
        var atten = _dust.Unred(Wave, 0.0, "calzetti");
        atten.Should().OnlyContain(v => Math.Abs(v - 1.0) < 1e-9);
    }
}
