using BAVCL;
using BAVCL.Core;
using BAVCL.Modules.Masking;
using BAVCL.Modules.Statistics;
using BAVCL.Modules.Structural;
using BAVCL.Types;
using Falcon.Core.Constants;
using Falcon.Core.Services.IO;
using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using System.Runtime.CompilerServices;

namespace Falcon.Core.Services.Instrument;

public class ResolutionService
{
    const float IsCloseAtol = 1e-8f;

    public ResolutionService(GPU gpu, Vector wave_in, Vector sres_in, bool log10_in = false)
    {
        this.gpu = gpu;

        if (wave_in.Length != sres_in.Length)
            throw new Exception("wave_in and sres_in must be of equal length");

        waveVec = wave_in;
        sresVec = sres_in;
        ReadOnlySpan<float> waveCpu = wave_in.RetrieveReadOnlySpan();
        ReadOnlySpan<float> sresCpu = sres_in.RetrieveReadOnlySpan();
        interpolator = new LinearInterpolation(waveCpu.ToArray(), sresCpu.ToArray());
        waveGpu = new Vector(gpu, interpolator.X, 0, false);
        sresGpu = new Vector(gpu, interpolator.Y, 0, false);

        log10 = log10_in;

        if (log10_in)
        {
            dv = Spectrum_velocity_scale(waveGpu);
            return;
        }

        dw = waveGpu.GetAt(1) - waveGpu.GetAt(0);
    }

    readonly GPU gpu;
    readonly Vector waveVec;
    readonly Vector sresVec;
    readonly Vector waveGpu;
    readonly Vector sresGpu;

    protected bool log10 { get; private set; }
    protected float dv { get; private set; }
    protected float dw { get; private set; }

    /// <summary>Velocity step per pixel in km/s when <see cref="log10"/> is true (firefly <c>dv</c>).</summary>
    public float Dv => dv;

    public LinearInterpolation interpolator { get; set; }
    protected float min_sig { get; set; }
    public float sig_vo { get; set; }
    public Vector sig_pd { get; set; }
    public Mask Sig_mask { get; set; }

    private float Spectrum_velocity_scale(Vector wave) =>
        FireflyConstants.c_kms * Spectral_coordinate_step(wave, log: true, _base: MathF.E);

    private static float Spectral_coordinate_step(Vector wave, bool log = false, float _base = 10f)
    {
        if (_base == MathF.E && log)
            return Vector.DiffX(wave.LogXIP(_base)).MeanX();

        if (log)
            return Vector.DiffX(wave).LogXIP(_base).IPOP(XMath.Rcp(XMath.Log(_base)), Operations.multiply).MeanX();

        return Vector.DiffX(wave).MeanX();
    }

    public void match(ResolutionService new_sres, bool no_offset = true, float min_sig_pix = 0f) =>
        GaussianKernelDifference(new_sres, no_offset, min_sig_pix);

    private void GaussianKernelDifference(ResolutionService new_sres, bool no_offset = true, float min_sig_pix = 0f)
    {
        min_sig = min_sig_pix;

        Vector _wave = waveGpu;
        Vector _sres = sresGpu;

        Vector interp_sres = new(gpu, new float[_wave.Length], 0, false);
        ReadOnlySpan<float> waveSpan = _wave.RetrieveReadOnlySpan();

        using (var scope = interp_sres.CpuScopeAndSync())
        {
            EditableView<float> view = scope.View;
            for (int i = 0; i < waveSpan.Length; i++)
                view[i] = new_sres.interpolator.Interpolate(waveSpan[i]);
        }

        Vector sig2_vd = DetermineVariance(_wave, interp_sres, _sres);

        if (no_offset)
        {
            sig_vo = 0f;

            if (log10)
            {
                finalize_GaussianKernelDifference(sig2_vd * XMath.Rcp(dv * dv));
                return;
            }

            finalize_GaussianKernelDifference(sig2_vd * XMath.Rcp(dw * dw));
            return;
        }

        float neg_amin_sig2_vd = -sig2_vd.MinX();
        float maxLocalDv = ComputeMaxLocalDv(_wave);
        sig_vo = neg_amin_sig2_vd - MathF.Pow(min_sig * maxLocalDv, 2f);

        if (sig_vo > 0f)
        {
            sig2_vd = sig2_vd + sig_vo;
            sig_vo = MathF.Sqrt(sig_vo);
        }
        else
        {
            sig_vo = 0f;
        }

        Vector sig2_pd = convert_vd2pd(sig2_vd, _wave);
        finalize_GaussianKernelDifference(sig2_pd);
    }

    private static float ComputeMaxLocalDv(Vector wave)
    {
        ReadOnlySpan<float> w = wave.RetrieveReadOnlySpan();
        float maxDv = 0f;
        for (int i = 0; i < w.Length - 1; i++)
        {
            float local = FireflyConstants.c_kms * (2f * (w[i + 1] - w[i]) / (w[i + 1] + w[i]));
            if (local > maxDv) maxDv = local;
        }
        return maxDv;
    }

    private static Vector DetermineVariance(Vector wave, Vector interp_sres, Vector sres)
    {
        Vector cOverWave = FireflyConstants.c_kms / wave;
        Vector cOverWaveSq = cOverWave ^ 2f;
        Vector waveOverFwhmSq = (wave / FireflyConstants.sig2FWHM) ^ 2f;
        Vector invInterpSq = interp_sres ^ (-2f);
        Vector invSresSq = sres ^ (-2f);
        return cOverWaveSq * waveOverFwhmSq * (invInterpSq - invSresSq);
    }

    private Vector convert_vd2pd(Vector sig2_vd, Vector wave)
    {
        if (log10)
            return sig2_vd * MathF.Pow(dv, -2f);

        float invCdwSq = XMath.Rcp(FireflyConstants.c_kms * FireflyConstants.c_kms * dw * dw);
        return (wave * wave) * sig2_vd * invCdwSq;
    }

    private void finalize_GaussianKernelDifference(Vector sig2_pd)
    {
        Mask zeroMask = sig2_pd.AbsX() <= IsCloseAtol;
        Vector absSig2 = sig2_pd.AbsX() ^ 0.5f;
        Vector normalized = sig2_pd / absSig2;
        sig_pd = normalized & (~zeroMask, 0f);
        Sig_mask = sig_pd < -min_sig;
    }

    public float[] adjusted_resolution(int[] indxs)
    {
        float sig2fwhmByCSq = MathF.Pow(FireflyConstants.sig2FWHM / FireflyConstants.c_kms, 2f);
        Vector sig2 = sig_pd * sig_pd.AbsX();

        if (indxs.Length == 0)
        {
            Vector pd2vd = convert_pd2vd(sig2);
            Vector inner = sig2fwhmByCSq * pd2vd + (sresGpu ^ (-2f));
            return ((inner ^ 0.5f).ReciprocalX()).RetrieveReadOnlySpan().ToArray();
        }

        float[] pd2vdArr = convert_pd2vd(sig2).RetrieveReadOnlySpan().ToArray();
        var output = new float[indxs.Length];
        for (int i = 0; i < indxs.Length; i++)
        {
            int idx = indxs[i];
            float selectedSres = 1f / MathF.Pow(interpolator.Y[idx], 2f);
            output[i] = 1f / MathF.Sqrt(sig2fwhmByCSq * MathF.Sqrt(pd2vdArr[idx]) + selectedSres);
        }
        return output;
    }

    private Vector convert_pd2vd(Vector sig2Pd)
    {
        if (log10)
            return sig2Pd * (dv * dv);

        float scale = FireflyConstants.c_kms * FireflyConstants.c_kms * dw * dw;
        return sig2Pd * scale / (waveGpu * waveGpu);
    }
}

class VariableGaussianKernel
{
    sealed class ResolutionGpuKernels(
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> multiplySumAxZero,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> kernCalc)
    {
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> MultiplySumAxZero { get; } = multiplySumAxZero;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int> KernCalc { get; } = kernCalc;
    }

    static readonly ConditionalWeakTable<GPU, ResolutionGpuKernels> KernelCache = new();

    static ResolutionGpuKernels GetKernels(GPU gpu) =>
        KernelCache.GetValue(gpu, static g =>
        {
            var acc = g.accelerator;
            return new ResolutionGpuKernels(
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int>(MultiplySumAxZeroKernel),
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int>(KernCalcKernel));
        });

    readonly ResolutionGpuKernels _kernels;
    public VariableGaussianKernel(GPU gpu, float[] sigma, float minsig = 0.01f, int nsig = 3, bool integral = false)
    {
        this.gpu = gpu;
        _kernels = GetKernels(gpu);
        n = sigma.Length;
        this.sigma = (float[])sigma.Clone();

        for (int i = 0; i < this.sigma.Length; i++)
            this.sigma[i] = MathF.Max(minsig, MathF.Abs(this.sigma[i]));

        p = (int)MathF.Ceiling(this.sigma.Max() * nsig);
        m = 2 * p + 1;

        float interval = (2f * MathF.Abs(p)) / m - 1;
        float[] x2 = Enumerable.Range(0, m).Select(val => MathF.Pow(-p + (val * interval), 2f)).ToArray();

        if (!integral)
            kernel = kernCalculation(gpu, new Vector(gpu, x2, 0, false), new Vector(gpu, this.sigma, 0, false));
        else
            Console.WriteLine("WARNING integral = true in Variable Gaussian KERNEL");
    }

    private int n { get; set; }
    private float[] sigma { get; set; }
    private int p { get; set; }
    private int m { get; set; }
    private Vector kernel { get; set; }
    readonly GPU gpu;

    public float[] Convolve(GPU gpu, float[] y) =>
        Convolve(gpu, new Vector(gpu, y, 0, false)).RetrieveReadOnlySpan().ToArray();

    public Vector Convolve(GPU gpu, Vector y)
    {
        if (y.Length != n)
            throw new Exception("Convolution ERROR : y array not the same length as n");

        using (GpuScope.BeginReadOnly(y, kernel))
        {
            Vector a = Create_a(y.RetrieveReadOnlySpan());
            return MultiplySumAxZero(a, kernel);
        }
    }

    Vector MultiplySumAxZero(Vector vectorA, Vector vectorB)
    {
        int columns = vectorA.Columns;
        int rows = vectorA.Length / columns;
        Vector output = new(gpu, new float[columns], 0, false);

        using (GpuScope.Begin(output, vectorA, vectorB))
        {
            _kernels.MultiplySumAxZero(
                gpu.accelerator.DefaultStream,
                columns,
                output.GetBuffer().View,
                vectorA.GetBuffer().View,
                vectorB.GetBuffer().View,
                columns,
                rows);
            gpu.accelerator.Synchronize();
        }

        return output;
    }

    static void MultiplySumAxZeroKernel(Index1D index, ArrayView<float> output, ArrayView<float> inputA, ArrayView<float> inputB, int columns, int rows)
    {
        float sum = 0f;
        for (int i = 0; i < rows; i++)
            sum += inputA[i * columns + index] * inputB[i * columns + index];
        output[index] = sum;
    }

    private Vector Create_a(ReadOnlySpan<float> y)
    {
        List<float> a = new();
        for (int i = 0; i < m; i++)
        {
            a.AddRange(new float[p]);
            for (int j = 0; j < n - m + 1; j++)
                a.Add(y[i + j]);
            a.AddRange(new float[p]);
        }
        return new Vector(gpu, a.ToArray(), kernel.Columns);
    }

    private Vector kernCalculation(GPU gpu, Vector x2, Vector sigmaVec)
    {
        int stride = sigmaVec.Columns > 1 ? sigmaVec.Columns : sigmaVec.Length;
        int kernelWidth = x2.Length;
        Vector output = new(gpu, new float[kernelWidth * stride], stride);

        using (GpuScope.Begin(output, x2, sigmaVec))
        {
            _kernels.KernCalc(
                gpu.accelerator.DefaultStream,
                stride,
                output.GetBuffer().View,
                x2.GetBuffer().View,
                sigmaVec.GetBuffer().View,
                kernelWidth,
                stride);
            gpu.accelerator.Synchronize();
        }

        return output;
    }

    static void KernCalcKernel(Index1D index, ArrayView<float> output, ArrayView<float> x2, ArrayView<float> sig, int kernelWidth, int stride)
    {
        float sum = 0f;
        for (int i = 0; i < kernelWidth; i++)
        {
            int indx = index + stride * i;
            output[indx] = XMath.Exp((-x2[i] * 0.5f) / XMath.Pow(sig[index], 2f));
            sum += output[indx];
        }

        sum = 1f / sum;
        for (int i = 0; i < kernelWidth; i++)
        {
            int indx = index + stride * i;
            output[indx] *= sum;
        }
    }
}
