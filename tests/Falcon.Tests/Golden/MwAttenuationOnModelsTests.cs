using BAVCL;
using Falcon.Core.Configuration;
using Falcon.Core.Constants;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.Dust;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Golden;

[Trait("Category", "Integration")]
public class MwAttenuationOnModelsTests
{
    [Fact]
    public void MwUnred_OnDowngradedScale_DoesNotZeroFlux()
    {
        string basePath = RepoRoot();
        string fits = Path.Combine(basePath, "assets", "sample-data", "spec-0266-51630-0016.fits");
        if (!File.Exists(fits)) return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(Path.Combine(basePath, "configs", "example.yaml"));
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, fits);
        var pipeline = new SpectralFitPipeline(gpu, fits, milky_Way_Reddening: true, hPF_Mode: false);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);

        var dust = new DustService(gpu);
        int nWave = FireflyConstants.wavelength.Length;
        var wave = new double[nWave];
        for (int k = 0; k < nWave; k++)
            wave[k] = FireflyConstants.wavelength.RetrieveReadOnlySpan()[k];

        double ebv = pipeline.Ebv_MW;
        var atten = dust.FitzpatrickUnred(wave, -ebv);
        double.IsFinite(atten[0]).Should().BeTrue("UV edge of MaStar grid must not yield infinite attenuation");

        int sidx = FireflyConstants.sidx;
        double raw = FireflyConstants.fluxgrid[10, 3, sidx, 100];
        double after = raw * atten[100];
        double.IsFinite(after).Should().BeTrue();
        after.Should().BeGreaterThan(0);
    }

    [Fact]
    public void MwCorrection_On0182_DoesNotDestroyModelFlux()
    {
        string basePath = RepoRoot();
        string fits = Path.Combine(basePath, "assets", "sample-data", "spec-0266-51630-0182.fits");
        if (!File.Exists(fits)) return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(Path.Combine(basePath, "configs", "example.yaml"));
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, fits);
        var pipeline = new SpectralFitPipeline(gpu, fits, milky_Way_Reddening: true, hPF_Mode: false);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);

        var dust = new DustService(gpu);
        int n = FireflyConstants.wavelength.Length;
        var wave = new double[n];
        ReadOnlySpan<float> w = FireflyConstants.wavelength.RetrieveReadOnlySpan();
        for (int i = 0; i < n; i++) wave[i] = w[i];

        double ebv = pipeline.Ebv_MW;
        var atten = dust.FitzpatrickUnred(wave, -ebv);
        atten.Min().Should().BeGreaterThan(0);
        atten.Should().OnlyContain(a => double.IsFinite(a) && a > 0);
    }

    static string RepoRoot() =>
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", ".."));
}
