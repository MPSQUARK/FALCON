using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.IO;
using BAVCL;
using Microsoft.Extensions.Logging;
using System.Diagnostics;

namespace Falcon.Core.Pipeline;

public sealed class SpectralFitRunner : ISpectralFitPipeline
{
    private readonly GPU _gpu;
    private readonly FalconFitSession _fitSession;
    private readonly IDustService _dust;
    private readonly IStatisticsService _statistics;
    private readonly ISpectralFitter _fitter;
    private readonly IFireflyFitsWriter _fitsWriter;
    private readonly ILogger<SpectralFitRunner> _logger;

    public SpectralFitRunner(
        GPU gpu,
        FalconFitSession fitSession,
        IDustService dust,
        IStatisticsService statistics,
        ISpectralFitter fitter,
        IFireflyFitsWriter fitsWriter,
        ILogger<SpectralFitRunner> logger)
    {
        _gpu = gpu;
        _fitSession = fitSession;
        _dust = dust;
        _statistics = statistics;
        _fitter = fitter;
        _fitsWriter = fitsWriter;
        _logger = logger;
    }

    public void FitFile(string fitsPath, CancellationToken cancellationToken = default)
    {
        if (_fitSession.IsShutdown)
            throw new InvalidOperationException("FALCON is shutting down.");

        FalconGpuHost.EnsureReady(_gpu);
        var fitToken = _fitSession.GetFitToken(cancellationToken);
        FitFileCore(fitsPath, fitToken);
    }

    public void FitBatch(IEnumerable<string> fitsPaths, CancellationToken cancellationToken = default, int maxParallelism = -1)
    {
        var paths = fitsPaths.ToArray();
        if (paths.Length == 0) throw new InvalidOperationException("No spectra to process.");

        int configuredMax = maxParallelism > 0 ? maxParallelism : LegacyConfig.MaxParallelFits;
        _logger.LogInformation(
            "Batch fit: {Count} spectra, parallelism {Parallelism}",
            paths.Length,
            FalconFitParallelism.DescribeParallelism(paths.Length, configuredMax));

        Parallel.ForEach(
            paths,
            FalconFitParallelism.BatchParallelOptions(paths.Length, configuredMax, cancellationToken),
            path => FitFile(path, cancellationToken));
    }

    private void FitFileCore(string fitsPath, CancellationToken cancellationToken)
    {
        string spectrum = Path.GetFileName(fitsPath);
        var sw = Stopwatch.StartNew();
        _logger.LogInformation("Starting fit for {Spectrum}", spectrum);
        try
        {
            cancellationToken.ThrowIfCancellationRequested();
            var (wavelength, flux, error, z, vdisp, ra, dec) = FitsSpectrumReader.ReadDataFits(_gpu, fitsPath);
            var pipeline = new SpectralFitPipeline(
                _gpu,
                fitsPath,
                LegacyConfig.Milky_Way_Reddening,
                LegacyConfig.HPF_Mode,
                LegacyConfig.N_Masked_Amstrongs,
                _dust, _statistics, _fitter, _fitsWriter,
                cancellationToken);
            pipeline.InitialiseSpectraParameters(wavelength, flux, error, z, [ra, dec], vdisp, LegacyConfig.Instrument_Resolution);
            pipeline.Fit_models_to_data();
            _logger.LogInformation("Finished fit for {Spectrum} in {ElapsedMs} ms", spectrum, sw.ElapsedMilliseconds);
        }
        catch (OperationCanceledException)
        {
            _logger.LogWarning("Fit cancelled for {Spectrum} after {ElapsedMs} ms", spectrum, sw.ElapsedMilliseconds);
            throw;
        }
        catch (Exception ex)
        {
            _logger.LogError(ex, "Fit failed for {Spectrum} after {ElapsedMs} ms", spectrum, sw.ElapsedMilliseconds);
            throw;
        }
    }
}
