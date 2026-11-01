using BAVCL;
using Falcon.Core.Pipeline;
using Microsoft.Extensions.Hosting;
using Microsoft.Extensions.Logging;

namespace Falcon.Core.DependencyInjection;

/// <summary>
/// Cancels in-flight fits and releases the ILGPU accelerator when the host stops.
/// </summary>
public sealed class FalconGpuLifetimeService(
    GPU gpu,
    FalconFitSession fitSession,
    ILogger<FalconGpuLifetimeService> logger) : IHostedService
{
    public Task StartAsync(CancellationToken cancellationToken) => Task.CompletedTask;

    public async Task StopAsync(CancellationToken cancellationToken)
    {
        logger.LogInformation("Shutting down FALCON GPU session");
        fitSession.RequestShutdown();

        // Give in-flight fits a brief window to observe cancellation before disposing the accelerator.
        await Task.Delay(TimeSpan.FromSeconds(2), cancellationToken).ConfigureAwait(false);

        try
        {
            gpu.Dispose();
            logger.LogInformation("GPU accelerator disposed");
        }
        catch (Exception ex)
        {
            logger.LogWarning(ex, "Error disposing GPU accelerator");
        }

        fitSession.Dispose();
    }
}
