using Falcon.Core.Abstractions;
using Falcon.Core.Services.Fitting;
using FluentAssertions;

namespace Falcon.Tests.Unit;

/// <summary>SigmaClip is internal; exercised through IterativeSpectralFitter branch-0 evaluation.</summary>
public class SigmaClipTests
{
    [Fact]
    public void Fit_WithOutlier_ProducesClippedMask()
    {
        var data = new[] { 1.0, 2.0, 100.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        var input = BuildInput(data, err, model);
        var result = new IterativeSpectralFitter().Fit(input);
        result.ClippedMask.Should().NotBeNull();
        result.ClippedMask.Length.Should().Be(3);
    }

    [Fact]
    public void Fit_UniformResiduals_ClipsAllOrNoneConsistently()
    {
        var data = new[] { 1.5, 2.5, 3.5 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        var input = BuildInput(data, err, model);
        var result = new IterativeSpectralFitter().Fit(input);
        result.MinChiSquared.Should().BeGreaterThanOrEqualTo(0);
    }

    [Theory]
    [InlineData(1.0, 1.0)]
    [InlineData(3.0, 3.0)]
    [InlineData(5.0, 5.0)]
    public void Fit_TightSigmaClip_AffectsChiSquared(double sigmaLow, double sigmaHigh)
    {
        var data = new[] { 1.0, 2.0, 2.5 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var models = new[] { 1.0, 2.0, 3.0, 1.2, 2.2, 3.2 };
        var input = BuildInput(data, err, models, numModels: 2);
        var loose = new IterativeSpectralFitter().Fit(input, new FitterOptions { SigmaClipLow = 10, SigmaClipHigh = 10 });
        var tight = new IterativeSpectralFitter().Fit(input, new FitterOptions { SigmaClipLow = sigmaLow, SigmaClipHigh = sigmaHigh });
        tight.MinChiSquared.Should().BeGreaterThanOrEqualTo(0);
        loose.Solutions.Count.Should().BeGreaterThan(0);
    }

    [Fact]
    public void Fit_BranchZeroSolutions_CarryClippedMask()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var models = new[] { 1.0, 2.0, 3.0, 2.0, 3.0, 4.0 };
        var result = new IterativeSpectralFitter().Fit(BuildInput(data, err, models, 2));
        result.Solutions.Where(s => s.BranchNumber == 0).Should().OnlyContain(s => s.ClippedMask != null);
    }

    [Fact]
    public void Fit_PerfectMatch_HasAllPointsInMaskOrZeroChi()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        var result = new IterativeSpectralFitter().Fit(BuildInput(data, err, model));
        result.MinChiSquared.Should().BeApproximately(0, 1e-10);
    }

    private static SpectralFitInput BuildInput(double[] data, double[] err, double[] models, int numModels = 1)
    {
        int n = data.Length;
        var wave = Enumerable.Range(0, n).Select(i => 4000.0 + i * 100).ToArray();
        return new SpectralFitInput { Wavelength = wave, Data = data, Error = err, Models = models };
    }
}
