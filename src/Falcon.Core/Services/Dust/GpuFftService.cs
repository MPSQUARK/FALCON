using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using BAVCL;
using BAVCL.Core;
using BAVCL.Core.Memory;
using Falcon.Core.Services.IO;
using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;

namespace Falcon.Core.Services.Dust;

/// <summary>Radix-2 FFT for HPF dust filtering. GPU kernels live in FALCON (not BAVCL).</summary>
public static class GpuFftService
{
    static readonly ConditionalWeakTable<GPU, FftKernels> KernelCache = new();
    static readonly ConditionalWeakTable<GPU, ConcurrentDictionary<(int N, int WindowSize, int WStart), HpfWindowVectors>> HpfWindowCache = new();

    public static double[] Hpf(GPU gpu, ReadOnlySpan<double> flux, int windowSize = 20, int wStart = 20)
    {
        ArgumentNullException.ThrowIfNull(gpu);

        int d = flux.Length;
        if (d == 0)
            return [];

        if (windowSize == 0 && wStart == 0)
            wStart = Math.Max(1, d / 100);

        int n = NextPowerOfTwo(d);
        wStart = Math.Clamp(wStart, 0, Math.Max(0, (n - 1) / 2));
        if (windowSize > 0)
            windowSize = Math.Min(windowSize, Math.Max(0, n - (2 * wStart)));

        var padded = PadToLength(flux, n);
        var windows = GetOrCreateHpfWindows(gpu, n, windowSize, wStart);
        return HpfGpu(gpu, padded, windows, d);
    }

    /// <summary>Forward FFT for unit tests. Length must be a power of two.</summary>
    internal static (float[] Real, float[] Imag) FftForward(GPU gpu, ReadOnlySpan<float> input)
    {
        ArgumentNullException.ThrowIfNull(gpu);
        int n = input.Length;
        if (n == 0 || (n & (n - 1)) != 0)
            throw new ArgumentException("FFT input length must be a positive power of two.", nameof(input));

        var src = new Vector(gpu, input.ToArray());
        using var scratch = RentFftScratch(gpu, n);
        var kernels = GetKernels(gpu);

        using (GpuScope.Begin(
            [scratch.RealA.PinTarget, scratch.ImagA.PinTarget, scratch.RealB.PinTarget, scratch.ImagB.PinTarget],
            [src]))
        {
            bool specInA = FftForwardGpu(
                gpu,
                kernels,
                src.GetBuffer().View,
                scratch.RealA.View,
                scratch.ImagA.View,
                scratch.RealB.View,
                scratch.ImagB.View,
                n,
                Log2(n));

            return specInA
                ? (PullReal(gpu, scratch.RealA), PullReal(gpu, scratch.ImagA))
                : (PullReal(gpu, scratch.RealB), PullReal(gpu, scratch.ImagB));
        }
    }

    /// <summary>Inverse FFT returning the real time-domain signal. Length must be a power of two.</summary>
    internal static float[] FftInverseReal(GPU gpu, ReadOnlySpan<float> real, ReadOnlySpan<float> imag)
    {
        ArgumentNullException.ThrowIfNull(gpu);
        int n = real.Length;
        if (n == 0 || imag.Length != n || (n & (n - 1)) != 0)
            throw new ArgumentException("FFT spectrum length must be a positive power of two.", nameof(real));

        var realVec = new Vector(gpu, real.ToArray());
        var imagVec = new Vector(gpu, imag.ToArray());
        using var scratch = RentFftScratch(gpu, n);
        var kernels = GetKernels(gpu);

        using (GpuScope.Begin(
            [scratch.RealA.PinTarget, scratch.ImagA.PinTarget, scratch.RealB.PinTarget, scratch.ImagB.PinTarget],
            [realVec, imagVec]))
        {
            CopySpectrumToA(gpu, kernels, realVec.GetBuffer().View, imagVec.GetBuffer().View, scratch, n);
            return FftInverseGpu(gpu, kernels, spectrumInA: true, scratch, n, Log2(n));
        }
    }

    static double[] HpfGpu(GPU gpu, double[] padded, HpfWindowVectors windows, int d)
    {
        int n = padded.Length;
        int log2n = Log2(n);
        var src = new Vector(gpu, ToFloat(padded));
        using var scratch = RentFftScratch(gpu, n);
        var kernels = GetKernels(gpu);

        using (GpuScope.Begin(
            [scratch.RealA.PinTarget, scratch.ImagA.PinTarget, scratch.RealB.PinTarget, scratch.ImagB.PinTarget],
            [src, windows.Window, windows.Unwindow]))
        {
            bool specInA = FftForwardGpu(
                gpu,
                kernels,
                src.GetBuffer().View,
                scratch.RealA.View,
                scratch.ImagA.View,
                scratch.RealB.View,
                scratch.ImagB.View,
                n,
                log2n);

            if (specInA)
            {
                FftComplexMulRealGpu(gpu, kernels, scratch.RealA.View, scratch.ImagA.View, windows.Window.GetBuffer().View, scratch.RealB.View, scratch.ImagB.View, n);
                var res = FftInverseGpu(gpu, kernels, spectrumInA: false, scratch, n, log2n);

                FftComplexMulRealGpu(gpu, kernels, scratch.RealB.View, scratch.ImagB.View, windows.Unwindow.GetBuffer().View, scratch.RealA.View, scratch.ImagA.View, n);
                var unres = FftInverseGpu(gpu, kernels, spectrumInA: false, scratch, n, log2n);

                return FinishHpf(res, unres, d);
            }

            FftComplexMulRealGpu(gpu, kernels, scratch.RealB.View, scratch.ImagB.View, windows.Window.GetBuffer().View, scratch.RealA.View, scratch.ImagA.View, n);
            var resB = FftInverseGpu(gpu, kernels, spectrumInA: true, scratch, n, log2n);

            FftComplexMulRealGpu(gpu, kernels, scratch.RealB.View, scratch.ImagB.View, windows.Unwindow.GetBuffer().View, scratch.RealA.View, scratch.ImagA.View, n);
            var unresB = FftInverseGpu(gpu, kernels, spectrumInA: true, scratch, n, log2n);

            return FinishHpf(resB, unresB, d);
        }
    }

    static double[] FinishHpf(float[] res, float[] unres, int d)
    {
        var resHead = new float[d];
        Array.Copy(res, resHead, d);
        float medRes = FitsSpectrumReader.Median(resHead);

        var output = new double[d];
        for (int i = 0; i < d; i++)
            output[i] = (1.0 + (res[i] - medRes) / (unres[i] == 0 ? 1e-30f : unres[i])) * medRes;
        return output;
    }

    static HpfWindowVectors GetOrCreateHpfWindows(GPU gpu, int n, int windowSize, int wStart)
    {
        var cache = HpfWindowCache.GetValue(gpu, static _ => new ConcurrentDictionary<(int, int, int), HpfWindowVectors>());
        return cache.GetOrAdd((n, windowSize, wStart), key =>
        {
            var window = BuildHpfWindow(key.N, key.WindowSize, key.WStart);
            var unwindow = new double[key.N];
            for (int i = 0; i < key.N; i++)
                unwindow[i] = 1 - window[i];
            unwindow[0] = 1;

            return new HpfWindowVectors(
                new Vector(gpu, ToFloat(window)),
                new Vector(gpu, ToFloat(unwindow)));
        });
    }

    static FftScratch RentFftScratch(GPU gpu, int n) =>
        new(
            BufferPools.For(gpu).Float.Rent(n),
            BufferPools.For(gpu).Float.Rent(n),
            BufferPools.For(gpu).Float.Rent(n),
            BufferPools.For(gpu).Float.Rent(n));

    static FftKernels GetKernels(GPU gpu) =>
        KernelCache.GetValue(gpu, static g =>
        {
            var accelerator = g.accelerator;
            return new FftKernels(
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(FftInitKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int>(FftBitReverseKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, float, int>(FftButterflyKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(FftScaleKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>>(NegateImagKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>>(FftComplexMulRealKernel),
                accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>>(CopyComplexKernel));
        });

    static void CopySpectrumToA(
        GPU gpu,
        FftKernels kernels,
        ArrayView<float> realSrc,
        ArrayView<float> imagSrc,
        FftScratch scratch,
        int n)
    {
        var stream = gpu.accelerator.DefaultStream;
        kernels.CopyComplex(stream, n, realSrc, imagSrc, scratch.RealA.View, scratch.ImagA.View);
        gpu.accelerator.Synchronize();
    }

    static bool FftForwardGpu(
        GPU gpu,
        FftKernels kernels,
        ArrayView<float> src,
        ArrayView<float> realA,
        ArrayView<float> imagA,
        ArrayView<float> realB,
        ArrayView<float> imagB,
        int n,
        int log2n)
    {
        var stream = gpu.accelerator.DefaultStream;
        kernels.Init(stream, n, src, realA, imagA);
        kernels.BitReverse(stream, n, realA, imagA, realB, imagB, log2n);
        RunButterflyStages(gpu, kernels, realB, imagB, realA, imagA, n, log2n, direction: 1f);
        gpu.accelerator.Synchronize();
        return (log2n % 2) != 0;
    }

    static float[] FftInverseGpu(
        GPU gpu,
        FftKernels kernels,
        bool spectrumInA,
        FftScratch scratch,
        int n,
        int log2n)
    {
        var stream = gpu.accelerator.DefaultStream;
        if (spectrumInA)
            kernels.NegateImag(stream, n, scratch.ImagA.View);
        else
            kernels.NegateImag(stream, n, scratch.ImagB.View);

        if (spectrumInA)
            kernels.BitReverse(stream, n, scratch.RealA.View, scratch.ImagA.View, scratch.RealB.View, scratch.ImagB.View, log2n);
        else
            kernels.BitReverse(stream, n, scratch.RealB.View, scratch.ImagB.View, scratch.RealA.View, scratch.ImagA.View, log2n);

        RunButterflyStages(gpu, kernels, scratch.RealB.View, scratch.ImagB.View, scratch.RealA.View, scratch.ImagA.View, n, log2n, direction: 1f);

        bool fftOutInA = (log2n % 2) != 0;
        float invN = 1f / n;
        if (fftOutInA)
        {
            kernels.NegateImag(stream, n, scratch.ImagA.View);
            kernels.Scale(stream, n, scratch.RealA.View, scratch.ImagA.View, invN);
            gpu.accelerator.Synchronize();
            return PullReal(gpu, scratch.RealA);
        }

        kernels.NegateImag(stream, n, scratch.ImagB.View);
        kernels.Scale(stream, n, scratch.RealB.View, scratch.ImagB.View, invN);
        gpu.accelerator.Synchronize();
        return PullReal(gpu, scratch.RealB);
    }

    static void RunButterflyStages(
        GPU gpu,
        FftKernels kernels,
        ArrayView<float> realA,
        ArrayView<float> imagA,
        ArrayView<float> realB,
        ArrayView<float> imagB,
        int n,
        int log2n,
        float direction)
    {
        var stream = gpu.accelerator.DefaultStream;
        int pairs = n >> 1;
        for (int stage = 0; stage < log2n; stage++)
        {
            bool toB = (stage % 2) == 0;
            if (toB)
                kernels.Butterfly(stream, pairs, realA, imagA, realB, imagB, stage, direction, n);
            else
                kernels.Butterfly(stream, pairs, realB, imagB, realA, imagA, stage, direction, n);
        }
    }

    static void FftComplexMulRealGpu(
        GPU gpu,
        FftKernels kernels,
        ArrayView<float> realIn,
        ArrayView<float> imagIn,
        ArrayView<float> window,
        ArrayView<float> realOut,
        ArrayView<float> imagOut,
        int n)
    {
        var stream = gpu.accelerator.DefaultStream;
        kernels.ComplexMulReal(stream, n, realIn, imagIn, window, realOut, imagOut);
        gpu.accelerator.Synchronize();
    }

    static float[] PullReal(GPU gpu, BufferEntity<float> entity)
    {
        using (GpuScope.Begin(entity.PinTarget))
        {
            gpu.accelerator.Synchronize();
            var host = new float[entity.Length];
            entity.View.CopyToCPU(host);
            return host;
        }
    }

    static double[] PadToLength(ReadOnlySpan<double> flux, int length)
    {
        var padded = new double[length];
        flux.CopyTo(padded);
        return padded;
    }

    static int NextPowerOfTwo(int n)
    {
        int pow2 = 1;
        while (pow2 < n)
            pow2 <<= 1;
        return pow2;
    }

    static int Log2(int n)
    {
        int log = 0;
        while ((1 << log) < n)
            log++;
        return log;
    }

    static float[] ToFloat(double[] values)
    {
        var output = new float[values.Length];
        for (int i = 0; i < values.Length; i++)
            output[i] = (float)values[i];
        return output;
    }

    static double[] BuildHpfWindow(int length, int windowSize, int wStart)
    {
        var window = new double[length];
        window[0] = 1;

        if (windowSize > 0)
        {
            double dw = windowSize;
            for (int i = 0; i < windowSize; i++)
            {
                int lo = wStart + i;
                int hi = length - 1 - (wStart + windowSize - i);
                if ((uint)lo < (uint)length)
                    window[lo] = (i + 1.0) / dw;
                if ((uint)hi < (uint)length)
                    window[hi] = (dw - i) / dw;
            }

            int passStart = wStart + windowSize;
            int passEnd = length - (wStart + windowSize);
            for (int i = passStart; i < passEnd; i++)
                window[i] = 1;
        }
        else
        {
            for (int i = wStart; i < length - wStart; i++)
                window[i] = 1;
        }

        return window;
    }

    readonly struct FftScratch(
        BufferEntity<float> realA,
        BufferEntity<float> imagA,
        BufferEntity<float> realB,
        BufferEntity<float> imagB) : IDisposable
    {
        public BufferEntity<float> RealA { get; } = realA;
        public BufferEntity<float> ImagA { get; } = imagA;
        public BufferEntity<float> RealB { get; } = realB;
        public BufferEntity<float> ImagB { get; } = imagB;

        public void Dispose()
        {
            RealA.Dispose();
            ImagA.Dispose();
            RealB.Dispose();
            ImagB.Dispose();
        }
    }

    sealed class FftKernels(
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>> init,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int> bitReverse,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, float, int> butterfly,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, float> scale,
        Action<AcceleratorStream, Index1D, ArrayView<float>> negateImag,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>> complexMulReal,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>> copyComplex)
    {
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>> Init { get; } = init;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int> BitReverse { get; } = bitReverse;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, float, int> Butterfly { get; } = butterfly;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, float> Scale { get; } = scale;
        public Action<AcceleratorStream, Index1D, ArrayView<float>> NegateImag { get; } = negateImag;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>> ComplexMulReal { get; } = complexMulReal;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>> CopyComplex { get; } = copyComplex;
    }

    sealed class HpfWindowVectors(Vector window, Vector unwindow)
    {
        public Vector Window { get; } = window;
        public Vector Unwindow { get; } = unwindow;
    }

    // KERNELS
    static void FftInitKernel(Index1D index, ArrayView<float> src, ArrayView<float> real, ArrayView<float> imag)
    {
        real[index] = src[index];
        imag[index] = 0f;
    }

    static void CopyComplexKernel(
        Index1D index,
        ArrayView<float> realSrc,
        ArrayView<float> imagSrc,
        ArrayView<float> realDst,
        ArrayView<float> imagDst)
    {
        realDst[index] = realSrc[index];
        imagDst[index] = imagSrc[index];
    }

    static void FftBitReverseKernel(
        Index1D index,
        ArrayView<float> realIn,
        ArrayView<float> imagIn,
        ArrayView<float> realOut,
        ArrayView<float> imagOut,
        int log2n)
    {
        int i = (int)index;
        int j = 0;
        int x = i;
        for (int b = 0; b < log2n; b++)
        {
            j = (j << 1) | (x & 1);
            x >>= 1;
        }

        realOut[j] = realIn[i];
        imagOut[j] = imagIn[i];
    }

    static void FftButterflyKernel(
        Index1D index,
        ArrayView<float> realIn,
        ArrayView<float> imagIn,
        ArrayView<float> realOut,
        ArrayView<float> imagOut,
        int stage,
        float direction,
        int n)
    {
        int half = 1 << stage;
        int step = half << 1;
        int pairIdx = (int)index;
        int group = pairIdx / half;
        int posInGroup = pairIdx % half;
        int i = group * step + posInGroup;
        int j = i + half;

        float theta = direction * (-2f * XMath.PI * posInGroup / step);
        float wr = XMath.Cos(theta);
        float wi = XMath.Sin(theta);

        float jr = realIn[j];
        float ji = imagIn[j];
        float tr = wr * jr - wi * ji;
        float ti = wr * ji + wi * jr;

        float ir = realIn[i];
        float ii = imagIn[i];
        realOut[i] = ir + tr;
        imagOut[i] = ii + ti;
        realOut[j] = ir - tr;
        imagOut[j] = ii - ti;
    }

    static void FftScaleKernel(Index1D index, ArrayView<float> real, ArrayView<float> imag, float scale)
    {
        real[index] *= scale;
        imag[index] *= scale;
    }

    static void NegateImagKernel(Index1D index, ArrayView<float> imag) =>
        imag[index] = -imag[index];

    static void FftComplexMulRealKernel(
        Index1D index,
        ArrayView<float> realIn,
        ArrayView<float> imagIn,
        ArrayView<float> window,
        ArrayView<float> realOut,
        ArrayView<float> imagOut)
    {
        float w = window[index];
        realOut[index] = realIn[index] * w;
        imagOut[index] = imagIn[index] * w;
    }
}
