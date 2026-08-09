using BAVCL;
using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.DependencyInjection;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Microsoft.Extensions.Logging;

var configPath = args.FirstOrDefault(a => a.EndsWith(".yaml", StringComparison.OrdinalIgnoreCase) || a.EndsWith(".yml", StringComparison.OrdinalIgnoreCase))
    ?? Path.Combine(AppContext.BaseDirectory, "configs", "example.yaml");

if (!File.Exists(configPath))
{
    Console.Error.WriteLine($"Config not found: {configPath}");
    return 1;
}

var options = YamlConfigurationExtensions.LoadFireflyYaml(configPath);
var basePath = AppContext.BaseDirectory;

var host = Host.CreateDefaultBuilder(args)
    .ConfigureServices(services => services.AddFalconCore())
    .Build();

var gpu = host.Services.GetRequiredService<GPU>();
var logger = host.Services.GetRequiredService<ILoggerFactory>().CreateLogger("Falcon.Cli");
ServiceCollectionExtensions.InitializeGpuData(gpu, options, basePath);

var pipeline = host.Services.GetRequiredService<ISpectralFitPipeline>();
var inputDir = Path.IsPathRooted(options.Input.SpectrumPath)
    ? options.Input.SpectrumPath
    : Path.Combine(basePath, options.Input.SpectrumPath);

string[] files = Directory.Exists(inputDir)
    ? Directory.GetFiles(inputDir, options.Input.InputGlob, SearchOption.AllDirectories)
    : File.Exists(inputDir) ? [inputDir] : [];

if (files.Length == 0)
{
    logger.LogError("No input spectra found at {Path}", inputDir);
    return 2;
}

var watch = System.Diagnostics.Stopwatch.StartNew();
pipeline.FitBatch(files);
watch.Stop();

logger.LogInformation("Galaxies analysed {Count}", files.Length);
logger.LogInformation("Time taken {Seconds:F2}s ({PerSpectrum:F2}s per spectrum)",
    watch.Elapsed.TotalSeconds,
    watch.Elapsed.TotalSeconds / files.Length);

return 0;
