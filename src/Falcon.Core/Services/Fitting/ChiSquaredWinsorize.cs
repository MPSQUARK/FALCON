namespace Falcon.Core.Services.Fitting;

using Falcon.Core.Numerics;

internal static class ChiSquaredWinsorize
{
    public static void InPlace(Span<float> work, double lowSigma, double highSigma, out float lower, out float upper)
    {
        lower = float.NegativeInfinity;
        upper = float.PositiveInfinity;
        if (work.IsEmpty) return;

        SimdSpanOps.ReplaceNonFinite(work, 1e6f);

        while (true)
        {
            var (mean, std) = SimdSpanOps.SampleMeanStd(work);
            lower = (float)(mean - lowSigma * std);
            upper = (float)(mean + highSigma * std);
            if (!SimdSpanOps.ClampInPlace(work, lower, upper))
                break;
        }
    }

    public static double ChiSquaredFromBlended(
        ReadOnlySpan<float> blended,
        ReadOnlySpan<bool> sharedMask,
        int branchNumber,
        Span<bool> clippedMask,
        double sigmaClipLow,
        double sigmaClipHigh)
    {
        if (branchNumber == 0)
        {
            Span<float> winsorized = stackalloc float[blended.Length > 512 ? 0 : blended.Length];
            float[]? heap = null;
            if (blended.Length > 512)
                heap = new float[blended.Length];
            Span<float> work = blended.Length > 512 ? heap!.AsSpan() : winsorized;
            blended.CopyTo(work);
            InPlace(work, sigmaClipLow, sigmaClipHigh, out float lower, out float upper);
            for (int w = 0; w < blended.Length; w++)
                clippedMask[w] = blended[w] > lower && blended[w] < upper;
            return SimdSpanOps.SumSquares(work);
        }

        double chiSq = SimdSpanOps.SumSquaresMasked(blended, sharedMask);
        return SimdSpanOps.CountMasked(sharedMask) == 0 ? double.NaN : chiSq;
    }
}
