using BAVCL;
using Falcon.Core.Configuration;
using Falcon.Core.Constants;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Golden;

[Trait("Category", "Integration")]
public class FluxGridLoadTests
{
    [Fact]
    public void MaStarFluxGrid_LoadsNonZeroValuesForValidModel()
    {
        string fitsModel = ResolveMaStarFits();
        if (!File.Exists(fitsModel)) return;

        string basePath = RepoRoot();

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(
            Path.GetFullPath(Path.Combine(basePath, "configs", "example.yaml")));
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        FireflyConstants.fluxgrid.Should().NotBeNull();
        int sidx = FireflyConstants.sidx;
        float raw = FireflyConstants.fluxgrid[10, 3, sidx, 100];
        raw.Should().BeGreaterThan(0, "FITS fluxgrid should retain MaStar flux values");

        float sanitized = raw <= -90f || !float.IsFinite(raw) ? 0f : raw;
        sanitized.Should().BeGreaterThan(0);
    }

    [Fact]
    public void SfdDustMaps_LoadNativeFitsImages()
    {
        string fitsModel = ResolveMaStarFits();
        if (!File.Exists(fitsModel)) return;

        string basePath = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..",
            "src", "Falcon.Cli", "bin", "Release", "net10.0"));
        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(
            Path.Combine(basePath, "configs", "example.yaml"));
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        FireflyConstants.ngp_dust.GetLength(0).Should().Be(4096);
        FireflyConstants.ngp_dust.GetLength(1).Should().Be(4096);
        FireflyConstants.sgp_dust.GetLength(0).Should().Be(4096);
        FireflyConstants.sgp_dust.GetLength(1).Should().Be(4096);
        FireflyConstants.ngp_dust.Cast<float>().Any(float.IsFinite).Should().BeTrue();
        FireflyConstants.sgp_dust.Cast<float>().Any(float.IsFinite).Should().BeTrue();
    }

    [Fact]
    public void Get_Model_WithDowngrade_NoMw_ProducesNonZeroBank()
    {
        string basePath = RepoRoot();
        string fits = Path.Combine(basePath, "assets", "SampleData", "spec-0266-51630-0016.fits");
        if (!File.Exists(fits)) return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(Path.Combine(basePath, "configs", "example.yaml"));
        opts.Dust.MilkyWayReddening = false;
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, fits);
        var pipeline = new SpectralFitPipeline(gpu, fits, milky_Way_Reddening: false, hPF_Mode: false);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);
        pipeline.Fit_models_to_data();

        string resultsPath = Path.Combine(basePath, "output", Path.GetFileNameWithoutExtension(fits) + ".results.json");
        if (!File.Exists(resultsPath)) return;
        var doc = Falcon.Core.Services.IO.FalconResultsReader.Load(resultsPath);
        doc.TopSolutions[0].ModelFlux.Max(Math.Abs).Should().BeGreaterThan(0);
    }

    [Fact]
    public void Get_Model_WithoutDowngrade_ProducesNonZeroBank()
    {
        string basePath = RepoRoot();
        string fits = Path.Combine(basePath, "assets", "SampleData", "spec-0266-51630-0016.fits");
        if (!File.Exists(fits) || !File.Exists(ResolveMaStarFits()))
            return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(Path.Combine(basePath, "configs", "example.yaml"));
        opts.Model.DowngradeModels = false;
        opts.Dust.MilkyWayReddening = false;
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, fits);
        var pipeline = new SpectralFitPipeline(gpu, fits, milky_Way_Reddening: false, hPF_Mode: false);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);
        pipeline.Fit_models_to_data();

        string resultsPath = Path.Combine(basePath, "output", Path.GetFileNameWithoutExtension(fits) + ".results.json");
        if (!File.Exists(resultsPath)) return;
        var doc = FalconResultsReader.Load(resultsPath);
        doc.TopSolutions[0].ModelFlux.Max(Math.Abs).Should().BeGreaterThan(0);
    }

    static string ResolveMaStarFits() =>
        Path.Combine(RepoRoot(), "Assets", "StellarPopulationModels", "MaStar_SSP_v0.2.fits.gz");

    static string RepoRoot() =>
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", ".."));
}
