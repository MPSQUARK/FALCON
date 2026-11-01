using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class UnredTests : IDisposable
{
    private readonly GPU _gpu = FalconGpuHost.Create();
    private readonly DustService _dust;

    public UnredTests() => _dust = new DustService(_gpu);

    public void Dispose() => _gpu.Dispose();

    [Fact]
    public void Unred_OutputLengthMatchesInputWavelengths()
    {
        var wave = new[] { 4000.0, 5000.0, 6000.0, 7000.0 };
        _dust.Unred(wave, 0.2).Length.Should().Be(wave.Length);
    }

    [Theory]
    [InlineData(0.05)]
    [InlineData(0.2)]
    [InlineData(0.8)]
    public void Unred_ValuesInPhysicalRange(double ebv)
    {
        var wave = Enumerable.Range(0, 20).Select(i => 3500.0 + i * 100).ToArray();
        var atten = _dust.Unred(wave, ebv);
        atten.Should().OnlyContain(v => v > 0 && v <= 1.0);
    }

    [Fact]
    public void Unred_UvWavelengthsMoreAttenuatedThanOptical()
    {
        var uv = _dust.Unred([3500.0], 0.3)[0];
        var optical = _dust.Unred([6500.0], 0.3)[0];
        uv.Should().BeLessThan(optical);
    }

    [Fact]
    public void Unred_FarRedApproachesUnity()
    {
        var atten = _dust.Unred([25000.0], 0.5)[0];
        atten.Should().BeApproximately(1.0, 0.01);
    }

    [Theory]
    [InlineData(4000)]
    [InlineData(5500)]
    [InlineData(9000)]
    public void Unred_MonotonicWithEbvAtFixedWave(double wave)
    {
        var a = _dust.Unred([wave], 0.1)[0];
        var b = _dust.Unred([wave], 0.4)[0];
        b.Should().BeLessThan(a);
    }
}
