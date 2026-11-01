using Falcon.Core.Services.Statistics;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class CalculateAveragesPdfTests
{
    private static readonly double[] Age = [1e9, 5e9, 10e9];
    private static readonly double[] Metal = [-1.0, 0.0, 0.5];
    private static readonly double[] Probs = [0.2, 0.5, 0.3];

    [Fact]
    public void CalculateAveragesPdf_ContainsAllExpectedKeys()
    {
        var (light, mass, unnorm) = BuildWeights();
        var av = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 200, 3e24, 1e-17);
        av.Keys.Should().Contain(new[]
        {
            "light_age", "light_metal", "mass_age", "mass_metal", "stellar_mass",
            "light_age_1_sig_plus", "light_age_1_sig_minus"
        });
    }

    [Fact]
    public void CalculateAveragesPdf_LightAgeNearWeightedMean()
    {
        var light = new[] { new[] { 1.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0 }, new[] { 0.0, 0.0, 1.0 } };
        var mass = light;
        var unnorm = light;
        var av = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 200, 3e24, 1e-17);
        av["light_age"].Should().BeApproximately(5.9e9, 2e9);
    }

    [Theory]
    [InlineData(100)]
    [InlineData(200)]
    public void CalculateAveragesPdf_ServiceMatchesStatisticsMath(int sampling)
    {
        var stats = new StatisticsService();
        var (light, mass, unnorm) = BuildWeights();
        var fromMath = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, sampling, 3e24, 1e-17);
        var fromService = stats.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 1e-17, 3e24, sampling);
        fromService["light_age"].Should().BeApproximately(fromMath["light_age"], 1e-3);
        fromService["stellar_mass"].Should().BeApproximately(fromMath["stellar_mass"], 1e-3);
    }

    [Fact]
    public void CalculateAveragesPdf_StellarMassScalesWithFluxUnits()
    {
        var (light, mass, unnorm) = BuildWeights();
        var low = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 50, 3e24, 1e-18);
        var high = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 50, 3e24, 1e-16);
        high["stellar_mass"].Should().BeApproximately(low["stellar_mass"] * 100, low["stellar_mass"] * 0.01);
    }

    [Fact]
    public void CalculateAveragesPdf_EachPropertyHasSixErrorKeys()
    {
        var (light, mass, unnorm) = BuildWeights();
        var av = StatisticsMath.CalculateAveragesPdf(Probs, light, mass, unnorm, Age, Metal, 100, 3e24, 1e-17);
        foreach (var key in new[] { "light_age", "light_metal", "mass_age", "mass_metal", "stellar_mass" })
        {
            av.Should().ContainKey($"{key}_1_sig_plus");
            av.Should().ContainKey($"{key}_3_sig_minus");
        }
    }

    private static (double[][] Light, double[][] Mass, double[][] Unnorm) BuildWeights()
    {
        var light = new[]
        {
            new[] { 0.6, 0.3, 0.1 },
            new[] { 0.2, 0.5, 0.3 },
            new[] { 0.1, 0.2, 0.7 }
        };
        var factors = new[] { 2.0, 1.5, 1.0 };
        var (unnorm, mass) = StatisticsMath.LightWeightsToMass(light, factors);
        return (light, mass, unnorm);
    }
}
