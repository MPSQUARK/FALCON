using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.Constants;
using Falcon.Core.Services.Dust;
using Falcon.Core.Services.IO;
using Microsoft.Extensions.Logging;

namespace Falcon.Core.Pipeline;

public sealed class SpectralFitRunner : ISpectralFitPipeline
{
    private readonly GPU _gpu;
    private readonly IDustService _dust;
    private readonly ILogger<SpectralFitRunner> _logger;

    public SpectralFitRunner(GPU gpu, IDustService dust, ILogger<SpectralFitRunner> logger)
    {
        _gpu = gpu;
        _dust = dust;
        _logger = logger;
    }

    public void FitFile(string fitsPath)
    {
        _logger.LogInformation("Fitting {Path}", fitsPath);
        var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(_gpu, fitsPath);
        var pipeline = new SpectralFitPipeline(fitsPath, LegacyConfig.Milky_Way_Reddening, LegacyConfig.HPF_Mode, LegacyConfig.N_Masked_Amstrongs, _dust);
        pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);
        pipeline.Fit_models_to_data();
    }

    public void FitBatch(IEnumerable<string> fitsPaths, int maxParallelism = -1)
    {
        var paths = fitsPaths.ToArray();
        if (paths.Length == 0) throw new InvalidOperationException("No spectra to process.");
        Parallel.ForEach(paths, new ParallelOptions { MaxDegreeOfParallelism = maxParallelism }, FitFile);
    }
}
