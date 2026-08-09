using Falcon.Core.Abstractions;
using Falcon.Core.Services.Fitting;
using FluentAssertions;

namespace Falcon.Tests.Component;

public class FitterTests
{
    [Fact]
    public void SingleModelUnitWeight_MatchesDirectChiSquared()
    {
        var data = new[] { 1.0, 2.0, 3.0 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var model = new[] { 1.0, 2.0, 3.0 };
        var input = new SpectralFitInput
        {
            Wavelength = new[] { 4000.0, 4100.0, 4200.0 },
            Data = data,
            Error = err,
            Models = model
        };

        var fitter = new IterativeSpectralFitter();
        var result = fitter.Fit(input);
        result.Solutions.Should().NotBeEmpty();
        result.MinChiSquared.Should().BeApproximately(0, 1e-10);
    }

    [Fact]
    public void TwoModelBlend_FindsImprovedSolution()
    {
        var data = new[] { 1.5, 2.5, 3.5 };
        var err = new[] { 0.1, 0.1, 0.1 };
        var m0 = new[] { 1.0, 2.0, 3.0 };
        var m1 = new[] { 2.0, 3.0, 4.0 };
        var models = m0.Concat(m1).ToArray();
        var input = new SpectralFitInput
        {
            Wavelength = new[] { 4000.0, 4100.0, 4200.0 },
            Data = data,
            Error = err,
            Models = models
        };

        var fitter = new IterativeSpectralFitter();
        var result = fitter.Fit(input);
        result.MinChiSquared.Should().BeLessThan(1.0);
    }
}
