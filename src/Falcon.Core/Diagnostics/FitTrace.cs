using System.Diagnostics;
using System.Text.Json;

namespace Falcon.Core.Diagnostics;

/// <summary>
/// Emits <c>FIT_TRACE\t{stage}\t{json}</c> lines when <c>FALCON_FIT_TRACE=1</c>.
/// Use <see cref="Scripts.compare_fit_traces"/> to diff against reference firefly output.
/// </summary>
public static class FitTrace
{
    static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = false,
        NumberHandling = System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals
    };

    public static bool Enabled =>
        string.Equals(Environment.GetEnvironmentVariable("FALCON_FIT_TRACE"), "1", StringComparison.Ordinal)
        || string.Equals(Environment.GetEnvironmentVariable("FALCON_FIT_TRACE"), "true", StringComparison.OrdinalIgnoreCase);

    public static void Emit(string stage, Dictionary<string, object> data)
    {
        if (!Enabled) return;
        Console.Out.WriteLine($"FIT_TRACE\t{stage}\t{JsonSerializer.Serialize(data, JsonOptions)}");
        Console.Out.Flush();
    }

    public static FitTraceSpan BeginSpan(string stage) => new(stage);

    public readonly struct FitTraceSpan : IDisposable
    {
        readonly string _stage;
        readonly Stopwatch _watch;

        internal FitTraceSpan(string stage)
        {
            _stage = stage;
            _watch = Stopwatch.StartNew();
        }

        public void Dispose()
        {
            if (!Enabled) return;
            Emit(_stage, new Dictionary<string, object> { ["elapsed_ms"] = _watch.Elapsed.TotalMilliseconds });
        }
    }

    public static Dictionary<string, object> Stats(string label, ReadOnlySpan<double> values)
    {
        if (values.IsEmpty)
            return new Dictionary<string, object> { [$"{label}_n"] = 0 };

        double min = double.PositiveInfinity;
        double max = double.NegativeInfinity;
        int finite = 0;
        int positive = 0;
        double sum = 0;
        foreach (double v in values)
        {
            if (!double.IsFinite(v)) continue;
            finite++;
            sum += v;
            if (v > 0) positive++;
            if (v < min) min = v;
            if (v > max) max = v;
        }

        double median = Median(values);
        return new Dictionary<string, object>
        {
            [$"{label}_n"] = values.Length,
            [$"{label}_finite"] = finite,
            [$"{label}_positive"] = positive,
            [$"{label}_min"] = finite > 0 ? min : 0,
            [$"{label}_max"] = finite > 0 ? max : 0,
            [$"{label}_median"] = median,
            [$"{label}_mean"] = finite > 0 ? sum / finite : 0,
        };
    }

    public static Dictionary<string, object> ModelGridStats(double[][] models, string prefix = "model0")
    {
        if (models.Length == 0 || models[0].Length == 0)
            return new Dictionary<string, object> { [$"{prefix}_n_models"] = models.Length };

        var stats = Stats(prefix, models[0]);
        stats[$"{prefix}_n_models"] = models.Length;
        return stats;
    }

    static double Median(ReadOnlySpan<double> data)
    {
        var finite = new List<double>(data.Length);
        foreach (double v in data)
            if (double.IsFinite(v)) finite.Add(v);
        if (finite.Count == 0) return 0;
        finite.Sort();
        int n = finite.Count;
        return n % 2 == 0 ? (finite[n / 2 - 1] + finite[n / 2]) / 2 : finite[n / 2];
    }
}
