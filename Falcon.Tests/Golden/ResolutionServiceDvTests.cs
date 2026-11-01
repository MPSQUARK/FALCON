using BAVCL;
using Falcon.Core.Configuration;
using Falcon.Core.Constants;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Instrument;
using FluentAssertions;
using Falcon.Core.Services.IO;

namespace Falcon.Tests.Golden;

[Trait("Category", "Integration")]
public class ResolutionServiceDvTests
{
    const double ExpectedDv = 69.02976392336477;
    const double Rtol = 1e-5;

    [Fact]
    public void Dv_OnMaStarWavelengthGrid_MatchesFireflyReference()
    {
        string fits = ResolveMaStarFits();
        if (!File.Exists(fits)) return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "configs", "example.yaml")));
        LegacyConfig.ApplyFrom(opts, RepoRoot());
        LegacyConfig.Setup(gpu, opts);

        FireflyConstants.wavelength.Should().NotBeNull();
        var sres = new Vector(gpu, Enumerable.Repeat(2000f, FireflyConstants.wavelength.Length).ToArray());

        var res = new ResolutionService(gpu, FireflyConstants.wavelength, sres, log10_in: true);
        ((double)res.Dv).Should().BeApproximately(ExpectedDv, ExpectedDv * Rtol);
    }

    [Fact]
    public void Dv_NumpyFormula_MatchesExpected()
    {
        string fits = ResolveMaStarFits();
        if (!File.Exists(fits)) return;

        double[] values = FitsAssetReader.ReadImage(fits, 2);
        var wave = values.Take(4563).Select(value => (float)value).ToArray();

        double step = 0;
        for (int i = 1; i < wave.Length; i++)
            step += Math.Log(wave[i]) - Math.Log(wave[i - 1]);
        step /= wave.Length - 1;
        double dv = FireflyConstants.c_kms * step;

        dv.Should().BeApproximately(ExpectedDv, ExpectedDv * Rtol);
    }

    static string ResolveMaStarFits() =>
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..",
            "Assets", "StellarPopulationModels", "MaStar_SSP_v0.2.fits.gz"));

    static string RepoRoot() =>
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", ".."));
}
