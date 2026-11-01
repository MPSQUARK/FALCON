using System.Text.Json;
using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;

// One-off generator: dotnet test --filter DumpGpuFftGolden
namespace Falcon.Tests.Unit;

public class GpuFftGoldenDump
{
    [Fact(Skip = "Run manually to regenerate golden_gpu_fft.json and golden_hpf.json")]
    public void DumpGpuFftGolden()
    {
        using var gpu = FalconGpuHost.Create();
        var mixed64 = BuildMixedSignal(64);
        var (real64, imag64) = GpuFftService.FftForward(gpu, mixed64);
        var round64 = GpuFftService.FftInverseReal(gpu, real64, imag64);

        var dc = Enumerable.Repeat(2.5f, 64).ToArray();
        var (dcReal, dcImag) = GpuFftService.FftForward(gpu, dc);

        var cosine = BuildCosine(1024, 11);
        var (cosReal, cosImag) = GpuFftService.FftForward(gpu, cosine);

        var rng = new Random(12345);
        var specReal = Enumerable.Range(0, 256).Select(_ => (float)(rng.NextDouble() * 2 - 1)).ToArray();
        var specImag = Enumerable.Range(0, 256).Select(_ => (float)(rng.NextDouble() * 2 - 1)).ToArray();
        var inverse256 = GpuFftService.FftInverseReal(gpu, specReal, specImag);

        var hpfFlux = Enumerable.Range(0, 64).Select(i => 1.0 + 0.05 * Math.Sin(i * 0.3)).ToArray();
        var hpf = new DustService(gpu).Hpf(hpfFlux);

        var payload = new Dictionary<string, object>
        {
            ["forward_mixed_64_real_8"] = ToDouble(real64, 8),
            ["forward_mixed_64_imag_8"] = ToDouble(imag64, 8),
            ["roundtrip_mixed_64_8"] = ToDouble(round64, 8),
            ["forward_dc_64_real_0"] = (double)dcReal[0],
            ["forward_dc_64_imag_0"] = (double)dcImag[0],
            ["forward_cosine_1024_mag_11"] = Magnitude(cosReal, cosImag, 11),
            ["forward_cosine_1024_mag_1013"] = Magnitude(cosReal, cosImag, 1013),
            ["inverse_random_256_8"] = ToDouble(inverse256, 8),
            ["hpf_64_first_8"] = hpf.Take(8).ToArray(),
            ["hpf_64_median"] = hpf.OrderBy(v => v).ElementAt(hpf.Length / 2)
        };

        foreach (int n in new[] { 64, 128, 1024, 4096 })
        {
            var input = BuildMixedSignal(n);
            var (real, imag) = GpuFftService.FftForward(gpu, input);
            var recovered = GpuFftService.FftInverseReal(gpu, real, imag);
            payload[$"roundtrip_mixed_{n}_max_err"] = MaxAbsError(input, recovered);
        }

        Console.WriteLine(JsonSerializer.Serialize(payload, new JsonSerializerOptions { WriteIndented = true }));
    }

    static float[] BuildMixedSignal(int n)
    {
        var output = new float[n];
        for (int i = 0; i < n; i++)
            output[i] = 0.7f
                + 0.2f * MathF.Sin(2f * MathF.PI * 3f * i / n)
                + 0.1f * MathF.Cos(2f * MathF.PI * 17f * i / n);
        return output;
    }

    static float[] BuildCosine(int n, int k)
    {
        var output = new float[n];
        for (int i = 0; i < n; i++)
            output[i] = MathF.Cos(2f * MathF.PI * k * i / n);
        return output;
    }

    static double[] ToDouble(float[] values, int count) =>
        values.Take(count).Select(v => (double)v).ToArray();

    static double MaxAbsError(float[] expected, float[] actual)
    {
        double max = 0;
        for (int i = 0; i < expected.Length; i++)
            max = Math.Max(max, Math.Abs(expected[i] - actual[i]));
        return max;
    }

    static double Magnitude(float[] real, float[] imag, int index) =>
        Math.Sqrt(real[index] * real[index] + imag[index] * imag[index]);
}
