using Falcon.Core.Abstractions;
using Falcon.Core.Services.Fitting;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class IterativeSpectralFitterExtendedTests
{
    [Fact]
    public void Fit_MatchesGoldenFixture()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_fitter.json");
        var input = new SpectralFitInput
        {
            Wavelength = new[] { 4000.0, 4100.0, 4200.0 },
            Data = new[] { 1.5, 2.5, 3.5 },
            Error = new[] { 0.1, 0.1, 0.1 },
            Models = new[] { 1.0, 2.0, 3.0, 2.0, 3.0, 4.0 }
        };
        var result = new IterativeSpectralFitter().Fit(input);
        result.MinChiSquared.Should().BeApproximately(golden.RootElement.GetProperty("min_chi_squared").GetDouble(), 1e-6);
        result.Converged.Should().Be(golden.RootElement.GetProperty("converged").GetBoolean());
        result.Solutions.Count.Should().Be(golden.RootElement.GetProperty("solution_count").GetInt32());
    }

    [Fact]
    public void Fit_InvalidInput_ThrowsOnValidate()
    {
        var input = new SpectralFitInput
        {
            Wavelength = new[] { 4000.0 },
            Data = new[] { 1.0, 2.0 },
            Error = new[] { 0.1, 0.1 },
            Models = new[] { 1.0, 2.0 }
        };
        var act = () => new IterativeSpectralFitter().Fit(input);
        act.Should().Throw<ArgumentException>();
    }

    [Fact]
    public void Fit_MaxSspComponents_CapsSpawnedBlends()
    {
        var result = new IterativeSpectralFitter().Fit(
            BuildBlendInput(numModels: 6),
            new FitterOptions { MaxIterations = 6, MaxSspComponents = 2, MixTopCount = 0 });

        foreach (var sol in result.Solutions)
            FitSolutionFingerprint.CountActiveComponents(sol.Weights).Should().BeLessThanOrEqualTo(2);
    }

    [Fact]
    public void Fit_MaxSspComponents_AllowsMultiComponentExploration()
    {
        const int n = 5;
        var data = Enumerable.Repeat(0.2, n).ToArray();
        var err = Enumerable.Repeat(0.01, n).ToArray();
        var wave = Enumerable.Range(0, n).Select(i => 4000.0 + i * 100).ToArray();
        var models = new double[n * n];
        for (int m = 0; m < n; m++)
        for (int w = 0; w < n; w++)
            models[m * n + w] = w == m ? 1.0 : 0.0;

        var input = new SpectralFitInput { Wavelength = wave, Data = data, Error = err, Models = models };
        var result = new IterativeSpectralFitter().Fit(
            input,
            new FitterOptions { MaxIterations = 4, MaxSspComponents = 5, MixTopCount = 0 });

        result.Solutions
            .Select(s => FitSolutionFingerprint.CountActiveComponents(s.Weights))
            .Max()
            .Should()
            .BeGreaterThanOrEqualTo(4);
    }

    [Fact]
    public void Fit_SolutionsAreUniqueRatioTuples()
    {
        var result = new IterativeSpectralFitter().Fit(BuildBlendInput(numModels: 3));
        var keys = result.Solutions
            .Select(s => FitSolutionFingerprint.FromWeights(s.Weights))
            .ToList();
        keys.Should().OnlyHaveUniqueItems();
    }

    [Fact]
    public void Fit_MaxIterationsZero_StillReturnsSolutions()
    {
        var input = BuildBlendInput();
        var result = new IterativeSpectralFitter().Fit(input, new FitterOptions { MaxIterations = 0 });
        result.Solutions.Should().NotBeEmpty();
    }

    [Theory]
    [InlineData(1)]
    [InlineData(5)]
    [InlineData(10)]
    public void Fit_MixTopCount_ControlsMixtureSolutions(int mixTop)
    {
        var input = BuildBlendInput();
        var result = new IterativeSpectralFitter().Fit(input, new FitterOptions { MixTopCount = mixTop });
        result.Solutions.Count.Should().BeGreaterThan(mixTop);
    }

    [Fact]
    public void Fit_SolutionsSortedByChiSquaredAfterMix()
    {
        var result = new IterativeSpectralFitter().Fit(BuildBlendInput());
        var chis = result.Solutions.Select(s => s.ChiSquared).ToArray();
        chis.Should().BeInAscendingOrder();
    }

    [Fact]
    public void Fit_ConvergesWithinDefaultIterations()
    {
        var result = new IterativeSpectralFitter().Fit(BuildBlendInput());
        result.IterationCount.Should().BeGreaterThan(0);
        result.IterationCount.Should().BeLessThanOrEqualTo(10);
    }

    [Theory]
    [InlineData(100)]
    [InlineData(500)]
    public void Fit_FitPerIterationCap_LimitsGrowth(int cap)
    {
        var input = BuildBlendInput(numModels: 4);
        var result = new IterativeSpectralFitter().Fit(input, new FitterOptions { FitPerIterationCap = cap, MaxIterations = 5 });
        result.Solutions.Count.Should().BeLessThanOrEqualTo(cap * 5 + 10);
    }

    [Fact]
    public void Fit_AllSolutionWeights_SumToUnity()
    {
        var result = new IterativeSpectralFitter().Fit(BuildBlendInput(numModels: 3));
        foreach (var sol in result.Solutions)
        {
            double sum = sol.Weights.Sum();
            sum.Should().BeApproximately(1.0, 1e-9);
        }
    }

    [Fact]
    public void Fit_SpawnTreeOnly_UsesDiscreteRationalWeights()
    {
        var result = new IterativeSpectralFitter().Fit(
            BuildBlendInput(numModels: 2),
            new FitterOptions { MaxIterations = 1, MixTopCount = 0 });

        var spawnOnly = result.Solutions.Where(s => s.BranchNumber > 0).ToList();
        spawnOnly.Should().NotBeEmpty();
        foreach (var sol in spawnOnly)
        {
            sol.Weights.Sum().Should().BeApproximately(1.0, 1e-9);
            foreach (double w in sol.Weights)
                IsSimpleRational(w).Should().BeTrue($"weight {w} should be a small rational fraction");
        }
    }

    static bool IsSimpleRational(double w, int maxDenom = 12)
    {
        if (w < 1e-9) return true;
        for (int q = 1; q <= maxDenom; q++)
        for (int p = 0; p <= q; p++)
            if (Math.Abs(w - (double)p / q) < 1e-6)
                return true;
        return false;
    }

    private static SpectralFitInput BuildBlendInput(int numModels = 2)
    {
        var data = new[] { 1.5, 2.5, 3.5, 4.5 };
        var err = Enumerable.Repeat(0.1, 4).ToArray();
        var wave = Enumerable.Range(0, 4).Select(i => 4000.0 + i * 100).ToArray();
        var models = Enumerable.Range(0, numModels)
            .SelectMany(m => Enumerable.Range(0, 4).Select(i => (double)(i + 1 + m)))
            .ToArray();
        return new SpectralFitInput { Wavelength = wave, Data = data, Error = err, Models = models };
    }
}
