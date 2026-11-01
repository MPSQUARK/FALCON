using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Instrument;
using FluentAssertions;

namespace Falcon.Tests.Golden;

public class VariableGaussianKernelTests
{
    [Fact]
    public void Convolve_PreservesPositiveFlux()
    {
        var gpu = FalconGpuHost.Create();
        int n = 32;
        float[] sigma = Enumerable.Range(0, n).Select(i => 1.5f + i * 0.01f).ToArray();
        var vgk = new VariableGaussianKernel(gpu, sigma);
        float[] flux = Enumerable.Range(1, n).Select(i => (float)i).ToArray();
        float[] outFlux = vgk.Convolve(gpu, flux);
        outFlux.Length.Should().Be(n);
        outFlux.Max().Should().BeGreaterThan(0);
    }

    [Fact]
    public void Convolve_HandlesNegativeSigmaWidths()
    {
        var gpu = FalconGpuHost.Create();
        int n = 32;
        float[] sigma = Enumerable.Range(0, n).Select(i => -1.5f - i * 0.01f).ToArray();
        var vgk = new VariableGaussianKernel(gpu, sigma);
        float[] flux = Enumerable.Range(1, n).Select(i => (float)i).ToArray();
        float[] outFlux = vgk.Convolve(gpu, flux);
        outFlux.Length.Should().Be(n);
        outFlux.Max().Should().BeGreaterThan(0);
    }

    [Fact]
    public void Convolve_HandlesFullMaStarLength()
    {
        var gpu = FalconGpuHost.Create();
        int n = 4563;
        var sigma = Enumerable.Repeat(2f, n).ToArray();
        var vgk = new VariableGaussianKernel(gpu, sigma);
        var flux = Enumerable.Range(0, n).Select(i => (float)(6e29 * (0.5 + 0.5 * Math.Sin(i * 0.01)))).ToArray();
        float[] outFlux = vgk.Convolve(gpu, flux);
        outFlux.Max(v => float.IsFinite(v) ? v : 0).Should().BeGreaterThan(0);
    }
}
