using BAVCL;
using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.DependencyInjection;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Serilog;
using Serilog.Events;

var configPath = args.FirstOrDefault(a => a.EndsWith(".yaml", StringComparison.OrdinalIgnoreCase) || a.EndsWith(".yml", StringComparison.OrdinalIgnoreCase))
    ?? Path.Combine(AppContext.BaseDirectory, "configs", "example.yaml");

if (!File.Exists(configPath))
{
    Console.Error.WriteLine($"Config not found: {configPath}");
    return 1;
}

Log.Logger = new LoggerConfiguration()
    .MinimumLevel.Information()
    .MinimumLevel.Override("Microsoft", LogEventLevel.Warning)
    .WriteTo.Console()
    .CreateLogger();

try
{
    var options = YamlConfigurationExtensions.LoadFireflyYaml(configPath);
    var basePath = AppContext.BaseDirectory;

    using var host = Host.CreateDefaultBuilder(args)
        .UseSerilog()
        .ConfigureServices(services => services.AddFalconCore())
        .Build();

    var gpu = host.Services.GetRequiredService<GPU>();
    ServiceCollectionExtensions.InitializeGpuData(gpu, options, basePath);

    var pipeline = host.Services.GetRequiredService<ISpectralFitPipeline>();
    var inputDir = Path.IsPathRooted(options.Input.SpectrumPath)
        ? options.Input.SpectrumPath
        : Path.Combine(basePath, options.Input.SpectrumPath);

    string[] files = Directory.Exists(inputDir)
        ? Directory.GetFiles(inputDir, options.Input.InputGlob, SearchOption.AllDirectories)
        : File.Exists(inputDir) ? [inputDir] : [];

    var spectrumArgs = args
        .Where(a => !a.EndsWith(".yaml", StringComparison.OrdinalIgnoreCase)
                 && !a.EndsWith(".yml", StringComparison.OrdinalIgnoreCase))
        .ToArray();
    if (spectrumArgs.Length > 0)
    {
        files = spectrumArgs
            .Select(p => Path.IsPathRooted(p) ? p : Path.Combine(basePath, p))
            .Where(File.Exists)
            .ToArray();
    }

    if (files.Length == 0)
    {
        Log.Error("No input spectra found at {InputPath}", inputDir);
        return 2;
    }

    pipeline.FitBatch(files);
    await host.StopAsync();
    return 0;
}
catch (Exception ex)
{
    Log.Fatal(ex, "FALCON CLI failed");
    return 1;
}
finally
{
    Log.CloseAndFlush();
}
