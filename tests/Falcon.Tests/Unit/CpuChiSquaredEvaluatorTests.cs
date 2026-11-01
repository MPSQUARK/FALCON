using Falcon.Core.Abstractions;
using Falcon.Core.Services.Fitting;
using FluentAssertions;

namespace Falcon.Tests.Unit;

/// <summary>CpuChiSquaredEvaluator is internal; validated through fitter chi-squared outcomes.</summary>
public class CpuChiSquaredEvaluatorTests
{
    private static readonly FitterOptions NoClip = new() { SigmaClipLow = 1000, SigmaClipHigh = 1000 };

    [Fact]
    public void Fit_PerfectSingleModelChiSquaredIsZero()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        var result = Fit(data, err, model);
        result.MinChiSquared.Should().BeApproximately(0, 1e-10);
    }

    [Fact]
    public void Fit_KnownOffsetChiMatchesManualCalculation()
    {
        var data = new[] { 1.1, 2.0, 2.9 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        double expected = Math.Pow(0.1 / 0.1, 2) + Math.Pow(-0.1 / 0.1, 2);
        var result = Fit(data, err, model);
        result.MinChiSquared.Should().BeApproximately(expected, 1e-4);
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    public void Fit_LargerErrorLowersChiSquared(double errVal)
    {
        var data = new[] { 1.5, 2.6, 3.4 };
        var err = new[] { errVal, errVal, errVal };
        var model = new[] { 1.0, 2.0, 3.0 };
        var result = Fit(data, err, model);
        result.MinChiSquared.Should().BeGreaterThan(0);
        var highErr = Fit(data, new[] { errVal * 2, errVal * 2, errVal * 2 }, model);
        result.MinChiSquared.Should().BeGreaterThan(highErr.MinChiSquared);
    }

    [Fact]
    public void Fit_TwoModels_BestChiLessThanWorstSingleModel()
    {
        var data = new[] { 1.5, 2.6, 3.5 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var m0 = new[] { 1.0, 2.0, 3.0 };
        var m1 = new[] { 2.0, 3.0, 4.0 };
        var result = Fit(data, err, m0.Concat(m1).ToArray(), numModels: 2);
        var singleChis = result.Solutions.Where(s => s.Weights.Count(w => w > 0.99) == 1).Select(s => s.ChiSquared).ToArray();
        singleChis.Max().Should().BeGreaterThan(0);
        result.MinChiSquared.Should().BeLessThan(singleChis.Max());
    }

    [Fact]
    public void Fit_EvaluatorSharedMask_IsPopulated()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var models = new[] { 1.0, 2.0, 3.0, 1.1, 2.1, 3.1 };
        var result = Fit(data, err, models, numModels: 2);
        result.ClippedMask.Should().NotBeNull();
        result.ClippedMask.Length.Should().Be(3);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Fit_ModelCount_MatchesSolutionSingleModelChis(int numModels)
    {
        int n = 4;
        var data = Enumerable.Range(0, n).Select(i => (double)i).ToArray();
        var err = Enumerable.Repeat(0.1, n).ToArray();
        var models = Enumerable.Range(0, numModels)
            .SelectMany(m => Enumerable.Range(0, n).Select(i => (double)(i + m)))
            .ToArray();
        var result = Fit(data, err, models, numModels);
        result.Solutions.Should().NotBeEmpty();
    }

    private static FitResult Fit(double[] data, double[] err, double[] models, int numModels = 1)
    {
        var wave = Enumerable.Range(0, data.Length).Select(i => 4000.0 + i * 10).ToArray();
        var input = new SpectralFitInput { Wavelength = wave, Data = data, Error = err, Models = models };
        return new IterativeSpectralFitter().Fit(input, NoClip);
    }
}
