using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class StatisticsIntegrateTests
{
    [Fact]
    public void Interpolate_ExactGridPoint_ReturnsFluxValue()
    {
        var wave = new[] { 4000.0, 4100.0, 4200.0 };
        var flux = new[] { 1.0, 2.0, 3.0 };
        StatisticsMath.Interpolate(wave, flux, 4100.0).Should().Be(2.0);
    }

    [Fact]
    public void Interpolate_Midpoint_ReturnsAverage()
    {
        var wave = new[] { 4000.0, 4200.0 };
        var flux = new[] { 1.0, 3.0 };
        StatisticsMath.Interpolate(wave, flux, 4100.0).Should().BeApproximately(2.0, 1e-9);
    }

    [Theory]
    [InlineData(3900)]
    [InlineData(4300)]
    public void Interpolate_OutOfRange_ClampsToEdge(double x)
    {
        var wave = new[] { 4000.0, 4100.0, 4200.0 };
        var flux = new[] { 1.0, 2.0, 3.0 };
        var v = StatisticsMath.Interpolate(wave, flux, x);
        v.Should().BeOneOf(1.0, 3.0);
    }

    [Fact]
    public void IntegrateToBins_ConstantFlux_ReturnsConstant()
    {
        var wave = new[] { 4000.0, 4010.0, 4020.0, 4030.0 };
        var flux = new[] { 2.0, 2.0, 2.0, 2.0 };
        var bisect = StatisticsMath.BisectArray(wave);
        var bins = StatisticsMath.IntegrateToBins(wave, flux, bisect);
        bins.Should().OnlyContain(v => Math.Abs(v - 2.0) < 0.01);
    }

    [Fact]
    public void IntegrateToBins_OutputLengthIsBisectMinusOne()
    {
        var wave = Enumerable.Range(0, 10).Select(i => 4000.0 + i * 10).ToArray();
        var flux = wave.Select(w => w * 1e-4).ToArray();
        var bisect = StatisticsMath.BisectArray(wave);
        StatisticsMath.IntegrateToBins(wave, flux, bisect).Length.Should().Be(bisect.Length - 1);
    }

    [Fact]
    public void IntegrateToBins_LinearSpectrum_ProducesMonotonicBins()
    {
        var wave = new[] { 0.0, 1.0, 2.0, 3.0 };
        var flux = new[] { 0.0, 1.0, 2.0, 3.0 };
        var bisect = new[] { 0.0, 0.5, 1.5, 2.5, 3.0 };
        var bins = StatisticsMath.IntegrateToBins(wave, flux, bisect);
        bins[^1].Should().BeGreaterThan(bins[0]);
        bins.Should().OnlyContain(v => !double.IsNaN(v));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(20)]
    public void IntegrateToBins_VariousGridSizes_ProducesFiniteValues(int n)
    {
        var wave = Enumerable.Range(0, n).Select(i => (double)i).ToArray();
        var flux = wave.Select(x => x * x).ToArray();
        var bisect = StatisticsMath.BisectArray(wave);
        var bins = StatisticsMath.IntegrateToBins(wave, flux, bisect);
        bins.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v));
    }
}
