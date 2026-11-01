using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class BisectArrayTests
{
    [Fact]
    public void BisectArray_MatchesGoldenFixture()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_bisect.json");
        var input = FireflyReferenceLoader.GetDoubleArray(golden, "input");
        var expected = FireflyReferenceLoader.GetDoubleArray(golden, "bisect");
        var actual = StatisticsMath.BisectArray(input);
        GoldenTestHelpers.AssertArraysEqual(actual, expected);
    }

    [Theory]
    [InlineData(new[] { 1.0, 2.0, 3.0 }, new[] { 1.0, 1.5, 2.5, 3.0 })]
    [InlineData(new[] { 10.0, 20.0 }, new[] { 10.0, 15.0, 20.0 })]
    [InlineData(new[] { 5.0 }, new[] { 5.0, 5.0 })]
    public void BisectArray_TheoryCases(double[] input, double[] expected)
    {
        var actual = StatisticsMath.BisectArray(input);
        GoldenTestHelpers.AssertArraysEqual(actual, expected);
    }

    [Fact]
    public void BisectArray_LengthIsInputPlusOne()
    {
        var input = new[] { 100.0, 200.0, 300.0, 400.0 };
        StatisticsMath.BisectArray(input).Length.Should().Be(input.Length + 1);
    }

    [Fact]
    public void BisectArray_EndpointsMatchInputExtrema()
    {
        var input = new[] { 4000.0, 4010.0, 4020.0 };
        var bisect = StatisticsMath.BisectArray(input);
        bisect[0].Should().Be(input[0]);
        bisect[^1].Should().Be(input[^1]);
    }

    [Fact]
    public void BisectArray_InteriorPointsAreMidpoints()
    {
        var input = new[] { 0.0, 10.0, 30.0 };
        var bisect = StatisticsMath.BisectArray(input);
        bisect[1].Should().Be(5.0);
        bisect[2].Should().Be(20.0);
    }
}
