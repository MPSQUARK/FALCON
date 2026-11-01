using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.Pipeline;
using Falcon.Core.Services.Cosmology;
using Falcon.Core.Services.Dust;
using Falcon.Core.Services.Fitting;
using Falcon.Core.Services.IO;
using Falcon.Core.Services.Instrument;
using Falcon.Core.Services.Statistics;
using Falcon.Core.Services.Manga;
using Falcon.Core.Services.Sfr;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using Microsoft.Extensions.Logging;

namespace Falcon.Core.DependencyInjection;

public static class ServiceCollectionExtensions
{
    public static IServiceCollection AddFalconCore(this IServiceCollection services)
    {
        services.AddSingleton<FalconFitSession>();
        services.AddSingleton<GPU>(_ => FalconGpuHost.Create());
        services.AddHostedService<FalconGpuLifetimeService>();
        services.AddSingleton<ICosmologyService, CosmologyServiceAdapter>();
        services.AddSingleton<ISpectrumReader, FitsSpectrumReaderAdapter>();
        services.AddSingleton<ISpectralFitter>(sp => new IterativeSpectralFitter(sp.GetRequiredService<GPU>()));
        services.AddSingleton<IDustService, DustService>();
        services.AddSingleton<IStatisticsService, StatisticsService>();
        services.AddSingleton<IFireflyFitsWriter, FireflyFitsWriter>();
        services.AddSingleton<IResolutionService, ResolutionServiceAdapter>();
        services.AddSingleton<ISpectralFitPipeline, SpectralFitRunner>();
        services.AddSingleton<IMangaBatchReader, MangaBatchReader>();
        services.AddSingleton<IStarFormationRateService, StarFormationRateService>();
        return services;
    }

    public static void InitializeGpuData(GPU gpu, FireflyOptions options, string basePath)
    {
        FalconGpuHost.EnsureReady(gpu);
        LegacyConfig.ApplyFrom(options, basePath);
        LegacyConfig.Setup(gpu, options);
    }
}

internal sealed class CosmologyServiceAdapter : ICosmologyService
{
    public float LuminosityDistance(GPU gpu, float redshift) => LegacyConfig.Cosmology.Luminosity_distance(gpu, redshift);
}

public interface ICosmologyService
{
    float LuminosityDistance(GPU gpu, float redshift);
}

internal sealed class FitsSpectrumReaderAdapter : ISpectrumReader
{
    public (Vector Wavelength, Vector Flux, Vector Error, float Redshift, float Vdisp, float Ra, float Dec) ReadFits(GPU gpu, string path) =>
        FitsSpectrumReader.ReadDataFits(gpu, path);
}

internal sealed class ResolutionServiceAdapter : IResolutionService;
