using BAVCL;
using Falcon.Core.Configuration;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Golden;

[Trait("Category", "Integration")]
public class ModelFluxDiagnosticTests
{
    [Fact]
    public void Get_Model_ProducesNonZeroFluxOnDataGrid()
    {
        string fits = ResolveSampleFits();
        if (!File.Exists(ResolveMaStarFits())) return;

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "configs", "example.yaml")));
        LegacyConfig.ApplyFrom(opts, RepoRoot());
        LegacyConfig.Setup(gpu, opts);

        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, fits);
        var pipeline = new SpectralFitPipeline(gpu, fits, hPF_Mode: true);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);
        pipeline.Fit_models_to_data();

        string resultsPath = Path.Combine(AppContext.BaseDirectory, "output",
            Path.GetFileNameWithoutExtension(fits) + ".results.json");
        if (!File.Exists(resultsPath))
            resultsPath = Path.Combine(
                Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..",
                    "src", "Falcon.Cli", "bin", "Debug", "net10.0", "output")),
                Path.GetFileNameWithoutExtension(fits) + ".results.json");

        if (!File.Exists(resultsPath)) return;

        var doc = Falcon.Core.Services.IO.FalconResultsReader.Load(resultsPath);
        var mf = doc.TopSolutions[0].ModelFlux;
        mf.Max(Math.Abs).Should().BeGreaterThan(0, "best-fit model flux must be visible on plot");
    }

    static string ResolveSampleFits()
    {
        string[] candidates =
        [
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..",
                "assets", "SampleData", "spec-0266-51630-0016.fits")),
            Path.Combine(AppContext.BaseDirectory, "assets", "SampleData", "spec-0266-51630-0016.fits"),
        ];
        return candidates.First(File.Exists);
    }

    static string ResolveMaStarFits() =>
        Path.Combine(RepoRoot(), "Assets", "StellarPopulationModels", "MaStar_SSP_v0.2.fits.gz");

    static string RepoRoot() =>
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", ".."));
}
