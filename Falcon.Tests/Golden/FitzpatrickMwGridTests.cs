using BAVCL;
using Falcon.Core.Configuration;
using Falcon.Core.Constants;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Golden;

[Trait("Category", "Integration")]
public class FitzpatrickMwGridTests
{
    [Fact]
    public void FitzpatrickUnred_OnMaStarGrid_IsFiniteAndPositive()
    {
        string basePath = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", ".."));

        var gpu = FalconGpuHost.Create();
        var opts = YamlConfigurationExtensions.LoadFireflyYaml(Path.Combine(basePath, "configs", "example.yaml"));
        LegacyConfig.ApplyFrom(opts, basePath);
        LegacyConfig.Setup(gpu, opts);

        int n = FireflyConstants.wavelength.Length;
        var wave = new double[n];
        ReadOnlySpan<float> w = FireflyConstants.wavelength.RetrieveReadOnlySpan();
        for (int i = 0; i < n; i++) wave[i] = w[i];

        var dust = new DustService(gpu);
        foreach (double ebv in new[] { 0.01, 0.05, 0.1, 0.5 })
        {
            var atten = dust.FitzpatrickUnred(wave, -ebv);
            atten.Should().AllSatisfy(a =>
            {
                double.IsFinite(a).Should().BeTrue($"ebv={ebv}");
                a.Should().BeGreaterThan(0, $"ebv={ebv}");
            });
            atten.Min().Should().BeGreaterThan(1e-6);
            atten.Max().Should().BeLessThan(1e6);
        }
    }
}
