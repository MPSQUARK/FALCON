using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class GpuFftTests : IDisposable
{
    private readonly GPU _gpu = FalconGpuHost.Create();

    public void Dispose() => _gpu.Dispose();

    [Fact]
    public void Forward_MixedSignal64_MatchesGoldenSpectrum()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        var input = BuildMixedSignal(64);
        var (real, imag) = GpuFftService.FftForward(_gpu, input);

        GoldenTestHelpers.AssertArraysEqual(
            real.Take(8).Select(v => (double)v).ToArray(),
            FireflyReferenceLoader.GetDoubleArray(golden, "forward_mixed_64_real_8"));
        GoldenTestHelpers.AssertArraysEqual(
            imag.Take(8).Select(v => (double)v).ToArray(),
            FireflyReferenceLoader.GetDoubleArray(golden, "forward_mixed_64_imag_8"));
    }

    [Fact]
    public void RoundTrip_MixedSignal64_MatchesGoldenSamples()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        var input = BuildMixedSignal(64);
        var (real, imag) = GpuFftService.FftForward(_gpu, input);
        var recovered = GpuFftService.FftInverseReal(_gpu, real, imag);

        GoldenTestHelpers.AssertArraysEqual(
            recovered.Take(8).Select(v => (double)v).ToArray(),
            FireflyReferenceLoader.GetDoubleArray(golden, "roundtrip_mixed_64_8"));
    }

    [Theory]
    [InlineData(64)]
    [InlineData(128)]
    [InlineData(1024)]
    [InlineData(4096)]
    public void RoundTrip_MixedSignal_MaxErrorBelowGolden(int n)
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        double expectedMaxErr = FireflyReferenceLoader.GetDouble(golden, $"roundtrip_mixed_{n}_max_err");

        var input = BuildMixedSignal(n);
        var (real, imag) = GpuFftService.FftForward(_gpu, input);
        var recovered = GpuFftService.FftInverseReal(_gpu, real, imag);

        double maxErr = 0;
        for (int i = 0; i < n; i++)
            maxErr = Math.Max(maxErr, Math.Abs(recovered[i] - input[i]));

        maxErr.Should().BeLessThanOrEqualTo(expectedMaxErr * 1.0001 + 1e-9);
    }

    [Fact]
    public void Forward_DcSignal_MatchesGoldenBinZero()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        var input = Enumerable.Repeat(2.5f, 64).ToArray();
        var (real, imag) = GpuFftService.FftForward(_gpu, input);

        real[0].Should().Be((float)FireflyReferenceLoader.GetDouble(golden, "forward_dc_64_real_0"));
        imag[0].Should().Be((float)FireflyReferenceLoader.GetDouble(golden, "forward_dc_64_imag_0"));

        for (int i = 1; i < 64; i++)
        {
            Math.Abs(real[i]).Should().BeLessThan(1e-5f);
            Math.Abs(imag[i]).Should().BeLessThan(1e-5f);
        }
    }

    [Fact]
    public void Forward_PureCosine_HasConjugateSymmetricPeaks()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        var input = BuildCosine(1024, 11);
        var (real, imag) = GpuFftService.FftForward(_gpu, input);

        Magnitude(real, imag, 11).Should().BeApproximately(
            FireflyReferenceLoader.GetDouble(golden, "forward_cosine_1024_mag_11"), 1e-3);
        Magnitude(real, imag, 1013).Should().BeApproximately(
            FireflyReferenceLoader.GetDouble(golden, "forward_cosine_1024_mag_1013"), 1e-3);
    }

    [Fact]
    public void Inverse_RandomSpectrum256_MatchesGoldenSamples()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_gpu_fft.json");
        var rng = new Random(12345);
        var real = Enumerable.Range(0, 256).Select(_ => (float)(rng.NextDouble() * 2 - 1)).ToArray();
        var imag = Enumerable.Range(0, 256).Select(_ => (float)(rng.NextDouble() * 2 - 1)).ToArray();

        var recovered = GpuFftService.FftInverseReal(_gpu, real, imag);
        GoldenTestHelpers.AssertArraysEqual(
            recovered.Take(8).Select(v => (double)v).ToArray(),
            FireflyReferenceLoader.GetDoubleArray(golden, "inverse_random_256_8"));
    }

    [Theory]
    [InlineData(64)]
    [InlineData(512)]
    [InlineData(4096)]
    public void Parseval_EnergyConserved(int n)
    {
        var input = BuildMixedSignal(n);
        var (real, imag) = GpuFftService.FftForward(_gpu, input);

        double timeEnergy = input.Sum(v => (double)v * v);
        double freqEnergy = 0;
        for (int i = 0; i < n; i++)
            freqEnergy += (real[i] * real[i] + imag[i] * imag[i]) / n;

        timeEnergy.Should().BeApproximately(freqEnergy, Math.Max(1e-4, 1e-6 * n));
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

    static double Magnitude(float[] real, float[] imag, int index) =>
        Math.Sqrt(real[index] * real[index] + imag[index] * imag[index]);
}
