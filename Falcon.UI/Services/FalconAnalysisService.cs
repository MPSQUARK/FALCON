using System.Collections.Concurrent;
using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.IO;
using BAVCL;

namespace Falcon.UI.Services;

public sealed class FalconAnalysisService
{
    private readonly ISpectralFitPipeline _pipeline;
    private readonly FalconFitSession _fitSession;
    private readonly GPU _gpu;
    private readonly IWebHostEnvironment _env;
    private readonly ILogger<FalconAnalysisService> _logger;

    public FalconAnalysisService(
        ISpectralFitPipeline pipeline,
        FalconFitSession fitSession,
        GPU gpu,
        IWebHostEnvironment env,
        ILogger<FalconAnalysisService> logger)
    {
        _pipeline = pipeline;
        _fitSession = fitSession;
        _gpu = gpu;
        _env = env;
        _logger = logger;
    }

    public string BasePath => AppContext.BaseDirectory;

    public string OutputDirectory => Path.Combine(BasePath, "output");

    public string PlotsDirectory
    {
        get
        {
            var dir = Path.Combine(_env.WebRootPath, "plots");
            Directory.CreateDirectory(dir);
            return dir;
        }
    }

    public GpuStatus GetGpuStatus()
    {
        try
        {
            FalconGpuHost.EnsureReady(_gpu);
            return new GpuStatus(true, _gpu.accelerator?.Name ?? "Accelerator ready");
        }
        catch (Exception ex)
        {
            return new GpuStatus(false, ex.Message);
        }
    }

    public FireflyOptions LoadConfig(string configRelativePath)
    {
        var path = ResolvePath(configRelativePath);
        return YamlConfigurationExtensions.LoadFireflyYaml(path);
    }

    public IReadOnlyList<SpectrumEntry> ListSpectra(string configRelativePath, IEnumerable<string>? extraPaths = null)
    {
        var opts = LoadConfig(configRelativePath);
        var paths = new HashSet<string>(DiscoverSpectra(opts), StringComparer.OrdinalIgnoreCase);
        if (extraPaths is not null)
        {
            foreach (var p in extraPaths.Where(File.Exists))
                paths.Add(Path.GetFullPath(p));
        }

        return paths
            .Select(path =>
            {
                string baseName = Path.GetFileNameWithoutExtension(path);
                string resultsPath = Path.Combine(OutputDirectory, baseName + ".results.json");
                return new SpectrumEntry(
                    path,
                    Path.GetFileName(path),
                    File.Exists(resultsPath),
                    resultsPath);
            })
            .OrderBy(s => s.FileName, StringComparer.OrdinalIgnoreCase)
            .ToList();
    }

    public async Task<RunSummary> RunBatchAsync(
        string configRelativePath,
        IReadOnlyList<string> selectedFiles,
        IProgress<BatchProgress>? progress = null,
        CancellationToken ct = default)
    {
        if (selectedFiles.Count == 0)
            throw new InvalidOperationException("Select at least one spectrum to analyse.");

        var opts = LoadConfig(configRelativePath);
        ServiceCollectionExtensions.InitializeGpuData(_gpu, opts, BasePath);

        var files = selectedFiles
            .Select(p => Path.IsPathRooted(p) ? p : Path.GetFullPath(ResolvePath(p)))
            .Where(File.Exists)
            .Distinct(StringComparer.OrdinalIgnoreCase)
            .ToArray();

        if (files.Length == 0)
            throw new InvalidOperationException("None of the selected spectra exist on disk.");

        if (files.Length < selectedFiles.Count)
            _logger.LogWarning("Skipped {Count} selected spectra that were not found", selectedFiles.Count - files.Length);

        var resultFiles = new ConcurrentBag<string>();
        var completedPaths = new ConcurrentBag<string>();
        using var linked = CancellationTokenSource.CreateLinkedTokenSource(_fitSession.ShutdownToken, ct);

        try
        {
            _logger.LogInformation(
                "UI batch: {Count} spectra, parallelism {Parallelism}",
                files.Length,
                FalconFitParallelism.DescribeParallelism(files.Length, opts.Fitting.MaxParallelFits));
            await Task.Run(() =>
            {
                Parallel.ForEach(
                    files,
                    FalconFitParallelism.BatchParallelOptions(files.Length, opts.Fitting.MaxParallelFits, linked.Token),
                    path =>
                    {
                        linked.Token.ThrowIfCancellationRequested();
                        progress?.Report(new BatchProgress(completedPaths.Count, files.Length, path, BatchPhase.Running));

                        try
                        {
                            _pipeline.FitFile(path, linked.Token);
                            string resultsPath = Path.Combine(OutputDirectory, Path.GetFileNameWithoutExtension(path) + ".results.json");
                            if (File.Exists(resultsPath))
                                resultFiles.Add(resultsPath);
                            completedPaths.Add(path);
                            progress?.Report(new BatchProgress(completedPaths.Count, files.Length, path, BatchPhase.Complete));
                        }
                        catch (OperationCanceledException)
                        {
                            throw;
                        }
                        catch (Exception ex)
                        {
                            _logger.LogError(ex, "Fit failed for {Spectrum}", Path.GetFileName(path));
                            completedPaths.Add(path);
                            progress?.Report(new BatchProgress(completedPaths.Count, files.Length, path, BatchPhase.Failed, ex.Message));
                        }
                    });
            }, linked.Token);
        }
        catch (OperationCanceledException)
        {
            throw;
        }

        return new RunSummary(files.Length, resultFiles.ToList());
    }

    public FalconSolutionRecord? LoadBestSolution(string resultsPath)
    {
        if (!File.Exists(resultsPath)) return null;
        try
        {
            var doc = FalconResultsReader.Load(resultsPath);
            return doc.TopSolutions.FirstOrDefault();
        }
        catch
        {
            return null;
        }
    }

    public FitResultSummary? LoadResultSummary(string resultsPath)
    {
        if (!File.Exists(resultsPath)) return null;
        try
        {
            var doc = FalconResultsReader.Load(resultsPath);
            var top = doc.TopSolutions.FirstOrDefault();
            doc.PdfAverages.TryGetValue("light_age", out double lightAge);
            doc.PdfAverages.TryGetValue("light_metal", out double lightMetal);
            double chi = top?.ChiSquared > 0 ? top.ChiSquared : doc.MinChiSquared;
            double chiDof = top?.ChiSquaredPerDof > 0
                ? top.ChiSquaredPerDof
                : (doc.DegreesOfFreedom > 0 ? doc.MinChiSquared / doc.DegreesOfFreedom : 0);
            double prob = top?.Probability ?? 0;
            return new FitResultSummary(
                chi,
                chiDof,
                prob,
                lightAge,
                lightMetal,
                top?.LightWeightedAgeGyr > 0 ? top.LightWeightedAgeGyr : lightAge,
                top?.LightWeightedMetal ?? 0);
        }
        catch
        {
            return null;
        }
    }

    private string[] DiscoverSpectra(FireflyOptions opts)
    {
        var inputDir = ResolvePath(opts.Input.SpectrumPath);
        return Directory.Exists(inputDir)
            ? Directory.GetFiles(inputDir, opts.Input.InputGlob, SearchOption.AllDirectories)
            : File.Exists(inputDir) ? [inputDir] : [];
    }

    public IReadOnlyList<string> ListResultFiles() =>
        FalconResultsReader.FindResultFiles(OutputDirectory).OrderByDescending(File.GetLastWriteTime).ToList();

    public FalconResultsDocument LoadResults(string resultsPath) =>
        FalconResultsReader.Load(resultsPath);

    public string GeneratePlot(string resultsJsonPath, int highlightRank = 1, PlotRangeMode rangeMode = PlotRangeMode.SnipToModelOverlap)
    {
        var doc = FalconResultsReader.Load(resultsJsonPath);
        string baseName = Path.GetFileNameWithoutExtension(resultsJsonPath);
        string plotName = $"{baseName}_r{highlightRank}.png";
        string plotPath = Path.Combine(PlotsDirectory, plotName);
        FalconPlotGenerator.SaveSpectrumPlot(doc, plotPath, highlightRank, rangeMode);
        return $"/plots/{plotName}";
    }

    public string ResolvePath(string path) =>
        Path.IsPathRooted(path) ? path : Path.Combine(BasePath, path);

    public static bool IsYamlPath(string path) =>
        path.EndsWith(".yaml", StringComparison.OrdinalIgnoreCase)
        || path.EndsWith(".yml", StringComparison.OrdinalIgnoreCase);

    public static bool IsSpectrumPath(string path) =>
        path.EndsWith(".fits", StringComparison.OrdinalIgnoreCase)
        || path.EndsWith(".dat", StringComparison.OrdinalIgnoreCase);
}

public sealed record GpuStatus(bool Ready, string Message);

public enum BatchPhase { Running, Complete, Failed }

public sealed record BatchProgress(int Completed, int Total, string CurrentFile, BatchPhase Phase, string? Error = null);

public sealed record RunSummary(int SpectraProcessed, IReadOnlyList<string> ResultFiles);

public sealed record SpectrumEntry(
    string FullPath,
    string FileName,
    bool HasExistingResults,
    string ResultsPath);

public sealed record FitResultSummary(
    double ChiSquared,
    double ChiSquaredPerDof,
    double Probability,
    double PdfLightAgeGyr,
    double PdfLightMetal,
    double TopLightAgeGyr,
    double TopLightMetal);
