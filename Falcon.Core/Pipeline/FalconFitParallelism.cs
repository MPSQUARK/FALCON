namespace Falcon.Core.Pipeline;

/// <summary>
/// Batch parallelism for multi-spectrum runs.
/// Matches <c>master</c> / <c>BAVCL_Update</c> <c>Program.cs</c>:
/// <c>Parallel.For(0, files.Length, i =&gt; {{ new Spectral_Model(...); Fit_models_to_data(); }})</c>
/// with no <c>MaxDegreeOfParallelism</c> cap (thread-pool default).
/// Set <see cref="LegacyConfig.MaxParallelFits"/> &gt; 0 only to limit concurrency.
/// </summary>
public static class FalconFitParallelism
{
    /// <summary>
    /// <c>MaxDegreeOfParallelism</c> for the batch loop.
    /// Returns -1 (unlimited) when <paramref name="configuredMax"/> is 0, matching legacy <c>Parallel.For</c>.
    /// </summary>
    public static int ResolveBatchParallelism(int fileCount, int configuredMax = 0)
    {
        if (fileCount <= 0) return 1;
        if (configuredMax > 0) return Math.Min(fileCount, configuredMax);
        return -1;
    }

    public static ParallelOptions BatchParallelOptions(int fileCount, int configuredMax, CancellationToken cancellationToken = default) =>
        new()
        {
            MaxDegreeOfParallelism = ResolveBatchParallelism(fileCount, configuredMax),
            CancellationToken = cancellationToken
        };

    public static string DescribeParallelism(int fileCount, int configuredMax)
    {
        int degree = ResolveBatchParallelism(fileCount, configuredMax);
        return degree < 0 ? $"unlimited (up to {fileCount} spectra)" : degree.ToString();
    }
}
