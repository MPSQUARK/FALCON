using System.Numerics;
using SimdF = System.Numerics.Vector<float>;
using SimdD = System.Numerics.Vector<double>;

namespace Falcon.Core.Numerics;

/// <summary>
/// Hardware-accelerated span math for CPU hot loops that are not on the GPU path.
/// Reductions use <see cref="double"/> with Kahan compensation so float32 inputs stay stable for large N.
/// </summary>
internal static class SimdSpanOps
{
    const float ChiFallback = 1e6f;
    const float MinError = 1e-30f;

    /// <summary>Compensated summation (Kahan). Use for all float32 → double reductions.</summary>
    static void KahanAdd(ref double sum, ref double compensation, double value)
    {
        double y = value - compensation;
        double t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
    }

    /// <summary>Float blend step with double intermediate so many model terms do not drift in f32.</summary>
    static void AddWeightedLane(ref float dest, float src, float weight) =>
        dest = (float)((double)dest + (double)weight * src);

    public static void AddWeighted(Span<float> dest, ReadOnlySpan<float> src, float weight)
    {
        if (weight == 0f) return;

        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            var w = new SimdF(weight);
            for (; i <= dest.Length - width; i += width)
            {
                var s = new SimdF(src.Slice(i, width));
                var products = s * w;
                for (int lane = 0; lane < width; lane++)
                    AddWeightedLane(ref dest[i + lane], products[lane], 1f);
            }
        }

        for (; i < dest.Length; i++)
            AddWeightedLane(ref dest[i], src[i], weight);
    }

    public static void AddInPlace(Span<float> dest, ReadOnlySpan<float> src) => AddWeighted(dest, src, 1f);

    public static void ScaleInPlace(Span<float> span, float scalar)
    {
        if (scalar == 1f || span.IsEmpty) return;

        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            var s = new SimdF(scalar);
            for (; i <= span.Length - width; i += width)
            {
                var v = new SimdF(span.Slice(i, width));
                (v * s).CopyTo(span.Slice(i, width));
            }
        }

        for (; i < span.Length; i++)
            span[i] *= scalar;
    }

    public static void ReplaceNonFinite(Span<float> span, float replacement)
    {
        for (int i = 0; i < span.Length; i++)
            if (!float.IsFinite(span[i])) span[i] = replacement;
    }

    public static double Sum(ReadOnlySpan<float> values)
    {
        double sum = 0;
        double compensation = 0;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            for (; i <= values.Length - width; i += width)
            {
                var v = new SimdF(values.Slice(i, width));
                for (int lane = 0; lane < width; lane++)
                    KahanAdd(ref sum, ref compensation, v[lane]);
            }
        }

        for (; i < values.Length; i++)
            KahanAdd(ref sum, ref compensation, values[i]);

        return sum;
    }

    public static (double Mean, double StdDev) SampleMeanStd(ReadOnlySpan<float> values)
    {
        if (values.IsEmpty) return (0, 0);
        double mean = Sum(values) / values.Length;

        double varSum = 0;
        double varCompensation = 0;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            for (; i <= values.Length - width; i += width)
            {
                var chunk = new SimdF(values.Slice(i, width));
                for (int lane = 0; lane < width; lane++)
                {
                    double d = chunk[lane] - mean;
                    KahanAdd(ref varSum, ref varCompensation, d * d);
                }
            }
        }

        for (; i < values.Length; i++)
        {
            double d = values[i] - mean;
            KahanAdd(ref varSum, ref varCompensation, d * d);
        }

        double std = values.Length > 1 ? Math.Sqrt(varSum / (values.Length - 1)) : 0;
        return (mean, std);
    }

    public static bool ClampInPlace(Span<float> span, float lower, float upper)
    {
        bool changed = false;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            var lo = new SimdF(lower);
            var hi = new SimdF(upper);
            for (; i <= span.Length - width; i += width)
            {
                var v = new SimdF(span.Slice(i, width));
                var clamped = global::System.Numerics.Vector.Min(global::System.Numerics.Vector.Max(v, lo), hi);
                if (!global::System.Numerics.Vector.EqualsAll(v, clamped))
                    changed = true;
                clamped.CopyTo(span.Slice(i, width));
            }
        }

        for (; i < span.Length; i++)
        {
            float v = span[i];
            if (v < lower) { span[i] = lower; changed = true; }
            else if (v > upper) { span[i] = upper; changed = true; }
        }

        return changed;
    }

    public static double SumSquares(ReadOnlySpan<float> span)
    {
        double sum = 0;
        double compensation = 0;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            for (; i <= span.Length - width; i += width)
            {
                var v = new SimdF(span.Slice(i, width));
                for (int lane = 0; lane < width; lane++)
                {
                    double x = v[lane];
                    KahanAdd(ref sum, ref compensation, x * x);
                }
            }
        }

        for (; i < span.Length; i++)
        {
            double x = span[i];
            KahanAdd(ref sum, ref compensation, x * x);
        }

        return sum;
    }

    public static double SumSquaresMasked(ReadOnlySpan<float> values, ReadOnlySpan<bool> mask)
    {
        double sum = 0;
        double compensation = 0;
        for (int i = 0; i < values.Length; i++)
        {
            if (!mask[i]) continue;
            double x = values[i];
            KahanAdd(ref sum, ref compensation, x * x);
        }

        return sum;
    }

    public static int CountMasked(ReadOnlySpan<bool> mask)
    {
        int count = 0;
        for (int i = 0; i < mask.Length; i++)
            if (mask[i]) count++;
        return count;
    }

    public static void ComputeChiRow(
        Span<float> dest,
        ReadOnlySpan<float> model,
        ReadOnlySpan<float> data,
        ReadOnlySpan<float> error)
    {
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            var fallback = new SimdF(ChiFallback);
            var minErr = new SimdF(MinError);
            for (; i <= dest.Length - width; i += width)
            {
                var m = new SimdF(model.Slice(i, width));
                var d = new SimdF(data.Slice(i, width));
                var e = global::System.Numerics.Vector.Max(new SimdF(error.Slice(i, width)), minErr);
                var chi = (m - d) / e;
                var bad = global::System.Numerics.Vector.GreaterThanOrEqual(global::System.Numerics.Vector.Abs(chi), new SimdF(float.MaxValue));
                global::System.Numerics.Vector.ConditionalSelect(bad, fallback, chi).CopyTo(dest.Slice(i, width));
            }
        }

        for (; i < dest.Length; i++)
        {
            float err = error[i];
            if (!float.IsFinite(err) || err <= 0f) err = MinError;
            float datum = float.IsFinite(data[i]) ? data[i] : 0f;
            float mod = float.IsFinite(model[i]) ? model[i] : 0f;
            float chi = (mod - datum) / err;
            dest[i] = float.IsFinite(chi) ? chi : ChiFallback;
        }
    }

    public static void PrepareErrorSpan(ReadOnlySpan<double> error, Span<float> errOut)
    {
        for (int i = 0; i < error.Length; i++)
        {
            double e = error[i];
            if (!double.IsFinite(e) || e <= 0) e = MinError;
            errOut[i] = (float)e;
        }
    }

    public static void PrepareDataSpan(ReadOnlySpan<double> data, Span<float> dataOut)
    {
        for (int i = 0; i < data.Length; i++)
        {
            double d = data[i];
            dataOut[i] = double.IsFinite(d) ? (float)d : 0f;
        }
    }

    public static void BuildChiBankRowMajor(
        Span<float> bank,
        ReadOnlySpan<float> data,
        ReadOnlySpan<float> error,
        ReadOnlySpan<double> models,
        Span<double> bankDouble,
        int modelIndex,
        int numWave)
    {
        int row = modelIndex * numWave;
        var modelF = new float[numWave];
        ReadOnlySpan<double> modelRow = models.Slice(row, numWave);
        for (int w = 0; w < numWave; w++)
            modelF[w] = double.IsFinite(modelRow[w]) ? (float)modelRow[w] : 0f;

        var chiRow = bank.Slice(row, numWave);
        ComputeChiRow(chiRow, modelF, data, error);
        for (int w = 0; w < numWave; w++)
            bankDouble[row + w] = chiRow[w];
    }

    public static void ScatterChiRowToGpuLayout(
        ReadOnlySpan<float> chiRow,
        Span<float> bank,
        int modelIndex,
        int numModels,
        int numWave)
    {
        for (int w = 0; w < numWave; w++)
            bank[w * numModels + modelIndex] = chiRow[w];
    }

    public static void BlendWeightedAccumulate(Span<double> dest, double weight, ReadOnlySpan<double> src)
    {
        if (weight == 0) return;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdD.Count;
            var w = new SimdD(weight);
            for (; i <= dest.Length - width; i += width)
            {
                var d = new SimdD(dest.Slice(i, width));
                var s = new SimdD(src.Slice(i, width));
                (d + s * w).CopyTo(dest.Slice(i, width));
            }
        }

        for (; i < dest.Length; i++)
            dest[i] += weight * src[i];
    }

    public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double sum = 0;
        double compensation = 0;
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdD.Count;
            for (; i <= a.Length - width; i += width)
            {
                var va = new SimdD(a.Slice(i, width));
                var vb = new SimdD(b.Slice(i, width));
                var products = va * vb;
                for (int lane = 0; lane < width; lane++)
                    KahanAdd(ref sum, ref compensation, products[lane]);
            }
        }

        for (; i < a.Length; i++)
            KahanAdd(ref sum, ref compensation, a[i] * b[i]);

        return sum;
    }

    public static void MultiplyFloatToDouble(
        ReadOnlySpan<float> flux,
        ReadOnlySpan<double> attenuation,
        Span<double> dest)
    {
        int i = 0;
        if (global::System.Numerics.Vector.IsHardwareAccelerated)
        {
            int width = SimdF.Count;
            for (; i <= flux.Length - width; i += width)
            {
                var f = new SimdF(flux.Slice(i, width));
                for (int lane = 0; lane < width; lane++)
                {
                    float v = f[lane];
                    if (!float.IsFinite(v) || v < 0f) v = 0f;
                    double scaled = (double)v * attenuation[i + lane];
                    dest[i + lane] = double.IsFinite(scaled) && scaled >= 0 ? scaled : 0;
                }
            }
        }

        for (; i < flux.Length; i++)
        {
            float v = flux[i];
            if (!float.IsFinite(v) || v < 0f) v = 0f;
            double scaled = (double)v * attenuation[i];
            dest[i] = double.IsFinite(scaled) && scaled >= 0 ? scaled : 0;
        }
    }

    public static float SanitizeFluxScalar(float value) =>
        !float.IsFinite(value) || value < 0f || value > float.MaxValue ? 0f : value;
}
