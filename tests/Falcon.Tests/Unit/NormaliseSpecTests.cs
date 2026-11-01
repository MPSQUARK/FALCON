using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class NormaliseSpecTests
{
    private readonly StatisticsService _stats = new();

    [Fact]
    public void NormaliseSpec_MatchesGoldenFixture()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_normalise.json");
        var data = new[] { 2.0, 2.0, 2.0, 2.0 };
        var models = new[] { new[] { 1.0, 1.0, 1.0, 1.0 }, new[] { 4.0, 4.0, 4.0, 4.0 } };
        var (normed, mass) = _stats.NormaliseSpec(data, models);
        GoldenTestHelpers.AssertArraysEqual(mass, FireflyReferenceLoader.GetDoubleArray(golden, "mass_factors"));
        GoldenTestHelpers.AssertArraysEqual(normed[0], FireflyReferenceLoader.GetDoubleArray(golden, "normed_model_0"));
        GoldenTestHelpers.AssertArraysEqual(normed[1], FireflyReferenceLoader.GetDoubleArray(golden, "normed_model_1"));
    }

    [Theory]
    [InlineData(2.0, 1.0, 2.0)]
    [InlineData(4.0, 2.0, 2.0)]
    [InlineData(1.0, 4.0, 0.25)]
    public void NormaliseSpec_MassFactorIsDataOverModelMedian(double dataVal, double modelVal, double expectedFactor)
    {
        var data = Enumerable.Repeat(dataVal, 4).ToArray();
        var models = new[] { Enumerable.Repeat(modelVal, 4).Select(v => (double)v).ToArray() };
        var (_, mass) = _stats.NormaliseSpec(data, models);
        mass[0].Should().BeApproximately(expectedFactor, 1e-9);
    }

    [Fact]
    public void NormaliseSpec_ScalesModelMedianToDataMedian()
    {
        var data = new[] { 3.0, 3.0, 3.0, 3.0 };
        var models = new[] { new[] { 1.0, 1.0, 1.0, 1.0 } };
        var (normed, _) = _stats.NormaliseSpec(data, models);
        normed[0].Average().Should().BeApproximately(3.0, 1e-9);
    }

    [Fact]
    public void NormaliseSpec_ZeroModelNormUsesIdentityMassFactor()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var models = new[] { new[] { 0.0, 0.0, 0.0 } };
        var (normed, mass) = _stats.NormaliseSpec(data, models);
        mass[0].Should().Be(1.0);
        normed[0].Should().Equal(models[0]);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(5)]
    public void NormaliseSpec_ReturnsOneMassFactorPerModel(int nModels)
    {
        var data = Enumerable.Repeat(2.0, 8).ToArray();
        var models = Enumerable.Range(0, nModels)
            .Select(m => Enumerable.Repeat(1.0 + m, 8).ToArray())
            .ToArray();
        var (normed, mass) = _stats.NormaliseSpec(data, models);
        normed.Length.Should().Be(nModels);
        mass.Length.Should().Be(nModels);
    }
}
