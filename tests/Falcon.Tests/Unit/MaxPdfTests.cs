using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class MaxPdfTests
{
    [Fact]
    public void MaxPdf_ConstantProperty_ReturnsUniformProbabilities()
    {
        var probs = new[] { 0.3, 0.5, 0.2 };
        var prop = new[] { 5.0, 5.0, 5.0 };
        var (xPdf, yPdf) = StatisticsMath.MaxPdf(probs, prop, sampling: 50);
        xPdf.Should().Equal(prop);
        yPdf.Should().OnlyContain(p => Math.Abs(p - 1.0 / 3.0) < 1e-9);
    }

    [Fact]
    public void MaxPdf_SpreadProperty_ReturnsNormalizedPdf()
    {
        var probs = new[] { 0.1, 0.9, 0.3 };
        var prop = new[] { 1.0, 2.0, 3.0 };
        var (_, yPdf) = StatisticsMath.MaxPdf(probs, prop, sampling: 100);
        yPdf.Sum().Should().BeApproximately(1.0, 1e-6);
        yPdf.Should().OnlyContain(v => v >= 0);
    }

    [Theory]
    [InlineData(10)]
    [InlineData(50)]
    [InlineData(200)]
    public void MaxPdf_SamplingParameter_ProducesValidLength(int sampling)
    {
        var probs = new[] { 0.5, 0.5 };
        var prop = new[] { 1.0, 2.0 };
        var (xPdf, yPdf) = StatisticsMath.MaxPdf(probs, prop, sampling);
        xPdf.Length.Should().BeGreaterThan(1);
        yPdf.Length.Should().Be(xPdf.Length - 2);
    }

    [Fact]
    public void MaxPdf_HighestProbBinDominatesNearPeakProperty()
    {
        var probs = new[] { 0.01, 0.98, 0.01 };
        var prop = new[] { 1.0, 5.0, 9.0 };
        var (xPdf, yPdf) = StatisticsMath.MaxPdf(probs, prop, sampling: 200);
        var peakIdx = yPdf.Select((v, i) => (v, i)).MaxBy(t => t.v).i;
        xPdf[peakIdx].Should().BeApproximately(5.0, 0.5);
    }

    [Fact]
    public void MaxPdf_PropertyPdfIsMonotonicInX()
    {
        var probs = new[] { 0.2, 0.6, 0.2 };
        var prop = new[] { 0.0, 1.0, 2.0 };
        var (xPdf, _) = StatisticsMath.MaxPdf(probs, prop, sampling: 100);
        for (int i = 1; i < xPdf.Length; i++)
            xPdf[i].Should().BeGreaterThanOrEqualTo(xPdf[i - 1]);
    }
}
