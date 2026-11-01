using System.Runtime.CompilerServices;
using BAVCL;
using ILGPU;
using ILGPU.Runtime;

namespace Falcon.Core.Services.Fitting;

internal static class ChiBlendGpu
{
    static readonly ConditionalWeakTable<GPU, ChiBlendKernels> KernelCache = new();

    internal sealed class ChiBlendKernels(
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> blendBatch)
    {
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> BlendBatch { get; } = blendBatch;
    }

    internal static ChiBlendKernels GetKernels(GPU gpu) =>
        KernelCache.GetValue(gpu, static g =>
        {
            var acc = g.accelerator;
            return new ChiBlendKernels(
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int>(BlendChiBatchKernel));
        });

    internal static void BlendBatch(
        GPU gpu,
        ArrayView<float> chiBank,
        ArrayView<float> weightRows,
        ArrayView<float> blendedOut,
        int batchSize,
        int numModels,
        int numWave)
    {
        var kernels = GetKernels(gpu);
        var stream = gpu.accelerator.DefaultStream;
        kernels.BlendBatch(stream, batchSize * numWave, blendedOut, chiBank, weightRows, numModels, numWave);
        gpu.accelerator.Synchronize();
    }

    static void BlendChiBatchKernel(
        Index1D index,
        ArrayView<float> blendedOut,
        ArrayView<float> chiBank,
        ArrayView<float> weightRows,
        int numModels,
        int numWave)
    {
        int b = index / numWave;
        int w = index % numWave;
        float sum = 0f;
        int weightRow = b * numModels;
        int chiRow = w * numModels;
        for (int m = 0; m < numModels; m++)
            sum += weightRows[weightRow + m] * chiBank[chiRow + m];
        blendedOut[b * numWave + w] = sum;
    }
}
