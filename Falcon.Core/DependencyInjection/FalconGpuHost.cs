using BAVCL;
using BAVCL.Core;

namespace Falcon.Core.DependencyInjection;

/// <summary>Creates and validates the application GPU used by FALCON pipelines.</summary>
public static class FalconGpuHost
{
    public static GPU Create()
    {
        var gpu = GPUManager.GetGPU();
        KernelModuleLoader.Load<float>(gpu, KernelWorkloads.Default);
        KernelModuleLoader.Load<float>(gpu, KernelWorkloads.Statistics);
        KernelModuleLoader.Load<float>(gpu, KernelWorkloads.Sorting);
        KernelModuleLoader.Load<int>(gpu, KernelWorkloads.Sorting);
        return gpu;
    }

    public static void EnsureReady(GPU gpu)
    {
        if (gpu is null)
            throw new InvalidOperationException("GPU instance is null.");

        if (gpu.accelerator is null)
            throw new InvalidOperationException(
                "GPU accelerator is not available. This often happens after a Blazor hot reload " +
                "left a disposed device in memory — restart the application and try again.");
    }
}
