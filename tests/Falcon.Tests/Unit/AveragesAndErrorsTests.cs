using Falcon.Core.Services.Statistics;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class AveragesAndErrorsTests
{
    [Fact]
    public void AveragesAndErrors_SingleProb_ReturnsSentinelValues()
    {
        var result = StatisticsMath.AveragesAndErrors([1.0], [5.0, 6.0], sampling: 50);
        result.Should().AllBeEquivalentTo(999.0);
    }

    [Fact]
    public void AveragesAndErrors_NaNProperty_ReturnsSentinelValues()
    {
        var probs = new[] { 0.5, 0.5 };
        var prop = new[] { double.NaN, 6.0 };
        var result = StatisticsMath.AveragesAndErrors(probs, prop, sampling: 50);
        result[0].Should().Be(999.0);
    }

    [Fact]
    public void AveragesAndErrors_ValidInput_ReturnsSevenStatistics()
    {
        var probs = new[] { 0.1, 0.7, 0.2 };
        var prop = new[] { 1.0, 5.0, 9.0 };
        var result = StatisticsMath.AveragesAndErrors(probs, prop, sampling: 200);
        result.Length.Should().Be(7);
        result[0].Should().NotBe(999.0);
    }

    [Theory]
    [InlineData(50)]
    [InlineData(100)]
    [InlineData(200)]
    public void AveragesAndErrors_UpperBoundsExceedLowerBounds(int sampling)
    {
        var probs = new[] { 0.2, 0.6, 0.2 };
        var prop = new[] { 2.0, 5.0, 8.0 };
        var stats = StatisticsMath.AveragesAndErrors(probs, prop, sampling);
        stats[1].Should().BeGreaterThanOrEqualTo(stats[0]);
        stats[2].Should().BeLessThanOrEqualTo(stats[0]);
        stats[3].Should().BeGreaterThanOrEqualTo(stats[1]);
        stats[4].Should().BeLessThanOrEqualTo(stats[2]);
    }

    [Fact]
    public void FindClosest_ReturnsNearestIndex()
    {
        var sorted = new[] { 0.0, 0.1, 0.5, 0.9, 1.0 };
        StatisticsMath.FindClosest(sorted, 0.48).Should().Be(2);
        StatisticsMath.FindClosest(sorted, 0.0).Should().Be(0);
    }

    [Fact]
    public void AveragesAndErrors_WiderSigmasSpanMoreValues()
    {
        var probs = new[] { 0.05, 0.9, 0.05 };
        var prop = new[] { 1.0, 10.0, 20.0 };
        var s = StatisticsMath.AveragesAndErrors(probs, prop, 200);
        (s[3] - s[4]).Should().BeGreaterThanOrEqualTo(s[1] - s[2]);
    }
}
