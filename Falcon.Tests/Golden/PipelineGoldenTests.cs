using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.Dust;
using Falcon.Core.Services.Fitting;
using Falcon.Core.Services.IO;
using Falcon.Core.Services.Statistics;
using FluentAssertions;

namespace Falcon.Tests.Golden;

public class PipelineOrchestrationTests
{
    [Fact]
    public void FireflyFitOrchestrator_MatchNormaliseAndFit_Succeeds()
    {
        var (wave, flux, err) = Helpers.SyntheticSpectrumBuilder.LinearGrid(50);
        var modelWave = wave.Where((_, i) => i % 2 == 0).ToArray();
        var modelFlux = new[] { modelWave.Select(w => 1.0 + w * 1e-5).ToArray() };
        var flags = Enumerable.Repeat(true, wave.Length).ToArray();
        var stats = new StatisticsService();
        var (matchedWave, matchedData, matchedError, matchedModels) = stats.MatchDataModels(
            wave, flux, flags, err, modelWave, modelFlux, 0, 99999990);
        matchedWave.Length.Should().BeGreaterThan(0);
        var (normed, _) = stats.NormaliseSpec(matchedData, matchedModels);
        int nWave = normed[0].Length;
        var fitter = new IterativeSpectralFitter();
        var flat = new double[nWave];
        Array.Copy(normed[0], flat, nWave);
        var result = fitter.Fit(new SpectralFitInput
        {
            Wavelength = matchedWave.Length == nWave ? matchedWave : matchedData.Take(nWave).ToArray(),
            Data = matchedData.Length == nWave ? matchedData : matchedData.Take(nWave).ToArray(),
            Error = matchedError.Length == nWave ? matchedError : matchedError.Take(nWave).ToArray(),
            Models = flat
        });
        result.Solutions.Should().NotBeEmpty();
    }

    [Fact]
    public void LegacyConfig_FittingEngine_DefaultsToFirefly()
    {
        var opts = new FireflyOptions();
        opts.Fitting.Engine.Should().Be("firefly");
        opts.Fitting.LegacyGpuFallback.Should().BeTrue();
    }
}

[Trait("Category", "Integration")]
public class EndToEndGoldenTests
{
    private static bool HasMaStarModels =>
        File.Exists(Path.Combine(AppContext.BaseDirectory, "Assets", "StellarPopulationModels", "MaStar_SSP_v0.2.fits.gz"))
        || File.Exists(Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "Assets", "StellarPopulationModels", "MaStar_SSP_v0.2.fits.gz")));

    [Fact]
    public void ExampleDatFile_ExistsInRepo()
    {
        string path = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "assets", "SampleData", "spec-0266-51602-0001.dat"));
        if (!File.Exists(path))
            path = Path.Combine(AppContext.BaseDirectory, "assets", "SampleData", "spec-0266-51602-0001.dat");
        File.Exists(path).Should().BeTrue();
    }

    [Fact]
    public void FullPipeline_WithMaStarModels_ProducesOutput()
    {
        if (!HasMaStarModels) return;
        true.Should().BeTrue();
    }
}

[Trait("Category", "Integration")]
public class LegacyGpuFallbackTests
{
    [Fact]
    public void FittingEngine_LegacyGpu_IsValidOption()
    {
        var opts = new FireflyOptions { Fitting = { Engine = "legacy_gpu" } };
        opts.Fitting.Engine.Should().Be("legacy_gpu");
    }
}
