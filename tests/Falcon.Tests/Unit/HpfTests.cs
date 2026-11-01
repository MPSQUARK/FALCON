using BAVCL;
using BAVCL.Modules.Statistics;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class HpfTests : IDisposable
{
    private readonly GPU _gpu = FalconGpuHost.Create();
    private readonly DustService _dust;

    public HpfTests() => _dust = new DustService(_gpu);

    public void Dispose() => _gpu.Dispose();

    [Fact]
    public void Hpf_MatchesGoldenFixtureFirstSamples()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_hpf.json");
        var flux = Enumerable.Range(0, 64).Select(i => 1.0 + 0.05 * Math.Sin(i * 0.3)).ToArray();
        var hpf = _dust.Hpf(flux);
        var expected = FireflyReferenceLoader.GetDoubleArray(golden, "hpf_first_8");
        GoldenTestHelpers.AssertArraysEqual(hpf.Take(8).ToArray(), expected);
    }

    [Fact]
    public void Hpf_PreservesLength()
    {
        var flux = Enumerable.Range(0, 128).Select(i => 1.0 + 0.01 * i).ToArray();
        _dust.Hpf(flux).Length.Should().Be(flux.Length);
    }

    [Fact]
    public void Hpf_MedianScaleApproximatelyPreserved()
    {
        var flux = Enumerable.Range(0, 128).Select(i => 1.0 + 0.01 * Math.Sin(i * 0.2)).ToArray();
        var hpf = _dust.Hpf(flux);
        new Vector(_gpu, hpf.Select(v => (float)v).ToArray()).MedianX()
            .Should().BeApproximately(new Vector(_gpu, flux.Select(v => (float)v).ToArray()).MedianX(), 0.5f);
    }

    [Theory]
    [InlineData(0, 0)]
    [InlineData(10, 10)]
    [InlineData(20, 20)]
    public void Hpf_WindowParameters_ProduceFiniteOutput(int windowSize, int wStart)
    {
        var flux = Enumerable.Repeat(1.0, 64).ToArray();
        var hpf = _dust.Hpf(flux, windowSize, wStart);
        hpf.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v));
    }

    [Fact]
    public void Hpf_ConstantInput_ReturnsNearConstant()
    {
        var flux = Enumerable.Repeat(3.14, 64).ToArray();
        var hpf = _dust.Hpf(flux);
        hpf.Should().OnlyContain(v => Math.Abs(v - 3.14) < 0.5);
    }

    [Fact]
    public void Hpf_ConstantInput_1024_ReturnsNearConstant()
    {
        var flux = Enumerable.Repeat(3.14, 1024).ToArray();
        var hpf = _dust.Hpf(flux);
        hpf.Should().OnlyContain(v => Math.Abs(v - 3.14) < 0.5);
    }

    [Theory]
    [InlineData(64)]
    [InlineData(128)]
    [InlineData(500)]
    public void Hpf_VariousLengths_ProduceFiniteOutput(int n)
    {
        var flux = Enumerable.Range(0, n).Select(i => 1.0 + 0.02 * Math.Cos(i * 0.1)).ToArray();
        _dust.Hpf(flux, windowSize: 10, wStart: 10)
            .Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v));
    }
}
