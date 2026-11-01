using Falcon.Core.Services.MassLoss;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class IsdEstimationTests
{
    private static readonly double[] X = [0.0, 1.0, 2.0, 3.0];
    private static readonly double[] Y = [0.0, 1.0, 0.0, 1.0];
    private static readonly double[] Z = [10.0, 20.0, 30.0, 40.0];

    [Fact]
    public void Estimate_AtKnownPoint_ReturnsExactValue()
    {
        var svc = new MassLossEstimationService(X, Y, Z);
        svc.Estimate(1.0, 1.0).Should().Be(20.0);
    }

    [Fact]
    public void Estimate_Midpoint_InterpolatesBetweenNeighbors()
    {
        var svc = new MassLossEstimationService(X, Y, Z);
        var value = svc.Estimate(0.5, 0.5);
        value.Should().BeInRange(10.0, 30.0);
    }

    [Theory]
    [InlineData(0.0, 0.0, 10.0)]
    [InlineData(2.0, 0.0, 30.0)]
    [InlineData(3.0, 1.0, 40.0)]
    public void Estimate_OnGridNodes(double x, double y, double expected)
    {
        new MassLossEstimationService(X, Y, Z).Estimate(x, y).Should().Be(expected);
    }

    [Fact]
    public void Estimate_OffGrid_WeightedByInverseDistanceSquared()
    {
        var svc = new MassLossEstimationService([0.0, 1.0], [0.0, 0.0], [0.0, 100.0]);
        var mid = svc.Estimate(0.5, 0.0);
        mid.Should().BeApproximately(50.0, 1.0);
    }

    [Fact]
    public void Estimate_FarFromPoints_StillFinite()
    {
        var svc = new MassLossEstimationService(X, Y, Z);
        svc.Estimate(10.0, 10.0).Should().BeGreaterThan(0);
    }

    [Theory]
    [InlineData(1.5, 0.5)]
    [InlineData(0.5, 1.5)]
    [InlineData(2.5, 0.5)]
    public void Estimate_VariousQueries_ReturnsWeightedAverage(double x, double y)
    {
        var svc = new MassLossEstimationService(X, Y, Z);
        var v = svc.Estimate(x, y);
        v.Should().BeInRange(10.0, 40.0);
    }
}
