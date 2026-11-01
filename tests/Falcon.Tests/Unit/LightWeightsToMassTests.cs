using Falcon.Core.Services.Statistics;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class LightWeightsToMassTests
{
    private readonly StatisticsService _stats = new();

    [Fact]
    public void LightWeightsToMassRow_NormalizesToSumOne()
    {
        var light = new[] { 0.25, 0.25, 0.5 };
        var factors = new[] { 2.0, 1.0, 0.5 };
        var mass = _stats.LightWeightsToMassRow(light, factors);
        mass.Sum().Should().BeApproximately(1.0, 1e-10);
    }

    [Fact]
    public void LightWeightsToMassRow_ZeroSum_ReturnsUnnormalized()
    {
        var light = new[] { 0.0, 0.0 };
        var factors = new[] { 1.0, 2.0 };
        var mass = StatisticsMath.LightWeightsToMassRow(light, factors);
        mass.Should().Equal(0.0, 0.0);
    }

    [Fact]
    public void LightWeightsToMass_MatrixVersion_NormalizesEachRow()
    {
        var light = new[] { new[] { 0.5, 0.5 }, new[] { 1.0, 0.0 } };
        var factors = new[] { 2.0, 1.0 };
        var (unnorm, mass) = _stats.LightWeightsToMass(light, factors);
        unnorm[0][0].Should().BeApproximately(1.0, 1e-9);
        mass[0].Sum().Should().BeApproximately(1.0, 1e-10);
        mass[1][0].Should().BeApproximately(1.0, 1e-9);
    }

    [Theory]
    [InlineData(1.0, 1.0, 1.0, 0.5)]
    [InlineData(0.5, 0.5, 2.0, 2.0 / 3.0)]
    public void LightWeightsToMassSingle_ReturnsExpectedMass(double w0, double w1, double f0, double expectedFirst)
    {
        var (sum, mass) = StatisticsMath.LightWeightsToMassSingle([w0, w1], [f0, 1.0]);
        sum.Should().BeApproximately(w0 * f0 + w1, 1e-9);
        mass[0].Should().BeApproximately(expectedFirst, 1e-9);
    }

    [Fact]
    public void LightWeightsToMass_UnnormMassMatchesWeightTimesFactor()
    {
        var light = new[] { new[] { 0.3, 0.7 } };
        var factors = new[] { 4.0, 2.0 };
        var (unnorm, _) = _stats.LightWeightsToMass(light, factors);
        unnorm[0][0].Should().BeApproximately(1.2, 1e-9);
        unnorm[0][1].Should().BeApproximately(1.4, 1e-9);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    public void LightWeightsToMass_PreservesRowCount(int nRows)
    {
        var light = Enumerable.Range(0, nRows).Select(_ => new[] { 0.5, 0.5 }).ToArray();
        var factors = new[] { 1.0, 1.0 };
        var (unnorm, mass) = _stats.LightWeightsToMass(light, factors);
        unnorm.Length.Should().Be(nRows);
        mass.Length.Should().Be(nRows);
    }
}
