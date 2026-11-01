namespace Falcon.Core.Services.Fitting;

/// <summary>
/// Canonical fingerprint for SSP weight vectors so near-duplicate mixes collapse to one solution.
/// </summary>
public static class FitSolutionFingerprint
{
    public const double MinFraction = 1e-4;
    public const int FractionDecimals = 4;

    public static string FromWeights(ReadOnlySpan<double> weights)
    {
        double sum = 0;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] > 0) sum += weights[i];
        if (sum <= 0) return "";

        var parts = new List<(int Index, double Frac)>();
        for (int i = 0; i < weights.Length; i++)
        {
            double frac = weights[i] / sum;
            if (frac >= MinFraction)
                parts.Add((i, frac));
        }

        parts.Sort((a, b) => a.Index.CompareTo(b.Index));
        return string.Join(",", parts.Select(p => $"{p.Index}:{Math.Round(p.Frac, FractionDecimals):F4}"));
    }

    public static int CountActiveComponents(ReadOnlySpan<double> weights)
    {
        double sum = 0;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] > 0) sum += weights[i];
        if (sum <= 0) return 0;

        int count = 0;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] / sum >= MinFraction)
                count++;
        return count;
    }
}
