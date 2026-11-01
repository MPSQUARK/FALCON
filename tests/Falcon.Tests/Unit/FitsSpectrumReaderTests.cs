using BAVCL;
using CFITS;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class FitsSpectrumReaderTests
{
    [Theory]
    [InlineData("spec-0266-51630-0016.fits", 0.12659f, 198.57f, 145.89217f, 0.059412f)]
    [InlineData("spec-0266-51630-0182.fits", 0.517111f, 425.94f, 145.89217f, 0.059412f)]
    public void ReadDataFits_MetadataMatchesFireflySdss(string fileName, float expectedZ, float expectedVdisp, float expectedRa, float expectedDec)
    {
        var gpu = FalconGpuHost.Create();
        var path = ResolveSampleFits(fileName);
        var (_, _, _, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, path);
        z.Should().BeApproximately(expectedZ, 0.001f);
        vdisp.Should().BeApproximately(expectedVdisp, 1f);
        ra.Should().BeApproximately(expectedRa, 0.001f);
        dec.Should().BeApproximately(expectedDec, 0.001f);
    }

    static string ResolveSampleFits(string fileName)
    {
        var candidates = new[]
        {
            Path.Combine(AppContext.BaseDirectory, "assets", "sample-data", fileName),
            Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "assets", "sample-data", fileName)),
        };
        return candidates.FirstOrDefault(File.Exists)
            ?? throw new FileNotFoundException($"Sample FITS not found: {fileName}");
    }

    static string SampleFitsPath() => ResolveSampleFits("spec-0266-51630-0016.fits");

    [Fact]
    public void SpecobjHdu_IsColumnTableWithNamedZColumn()
    {
        using var document = FitsDocument.Open(SampleFitsPath());
        FitsTable? table = document.Hdus[2].Table;
        table.Should().NotBeNull();
        FitsColumn z = table.Column("Z");
        z.ReadReal(0, 0).Should().BeApproximately(0.12659, 0.001);
    }

    [Fact]
    public void RawFitsTable_HasPositiveFluxPixels()
    {
        using var document = FitsDocument.Open(SampleFitsPath());
        FitsTable? table = document.Hdus[1].Table;
        table.Should().NotBeNull();
        FitsColumn fluxColumn = table.Column("flux");
        FitsColumn ivarColumn = table.Column("ivar");
        var flux = Enumerable.Range(0, table.RowCount)
            .Select(row => (float)fluxColumn.ReadReal(row, 0))
            .ToArray();
        var ivar = Enumerable.Range(0, table.RowCount)
            .Select(row => (float)ivarColumn.ReadReal(row, 0))
            .ToArray();

        flux.Count(x => x > 0).Should().BeGreaterThan(100, "sample SDSS flux should have many positive pixels");
        ivar.Count(x => x > 0).Should().BeGreaterThan(100);
    }

    [Fact]
    public void ReadDataFits_LeavesGoodPixelsAfterMasking()
    {
        var gpu = FalconGpuHost.Create();
        var (wave, flux, error, _, _, _, _) = FitsSpectrumReader.ReadDataFits(gpu, SampleFitsPath());

        var fluxCpu = flux.Pull();
        var errCpu = error.Pull();
        int good = 0;
        for (int i = 0; i < fluxCpu.Length; i++)
        {
            bool ok = !(float.IsNaN(fluxCpu[i]) || float.IsInfinity(fluxCpu[i]) || fluxCpu[i] <= 0
                || float.IsNaN(errCpu[i]) || float.IsInfinity(errCpu[i]));
            if (ok) good++;
        }

        good.Should().BeGreaterThan(100, $"got flux range [{fluxCpu.Min():G4}, {fluxCpu.Max():G4}]");
    }

    [Fact]
    public void InitialiseSpectraParameters_DoesNotRejectEntireSpectrum()
    {
        var gpu = FalconGpuHost.Create();
        var path = SampleFitsPath();
        var (wave, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(gpu, path);
        var pipeline = new SpectralFitPipeline(
            gpu,
            path,
            milky_Way_Reddening: false,
            hPF_Mode: false,
            n_Masked_Amstrongs: 20);

        var act = () => pipeline.InitialiseSpectraParameters(wave, flux, error, z, [ra, dec], vdisp, 2.5f);
        act.Should().NotThrow();
        pipeline.Flux.Length.Should().BeGreaterThan(100);
    }
}
