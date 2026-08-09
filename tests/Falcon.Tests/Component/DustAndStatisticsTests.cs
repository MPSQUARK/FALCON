using Falcon.Core.Services.Dust;
using Falcon.Core.Services.Statistics;
using FluentAssertions;

namespace Falcon.Tests.Component;

public class DustAndStatisticsTests
{
    [Fact]
    public void Hpf_PreservesMedianScale()
    {
        var flux = Enumerable.Range(0, 128).Select(i => 1.0 + 0.01 * Math.Sin(i * 0.2)).ToArray();
        var dust = new DustService();
        var hpf = dust.Hpf(flux);
        hpf.Length.Should().Be(flux.Length);
        hpf.Average().Should().BeApproximately(flux.Average(), 0.5);
    }

    [Fact]
    public void NormaliseSpec_ScalesModelsToDataMedian()
    {
        var data = new[] { 2.0, 2.0, 2.0, 2.0 };
        var models = new[] { new[] { 1.0, 1.0, 1.0, 1.0 }, new[] { 4.0, 4.0, 4.0, 4.0 } };
        var stats = new StatisticsService();
        var (normed, mass) = stats.NormaliseSpec(data, models);
        normed[0].Average().Should().BeApproximately(2.0, 1e-6);
        mass[0].Should().BeApproximately(2.0, 1e-6);
    }

    [Fact]
    public void ConvertChisToProbs_SumsToOne()
    {
        var stats = new StatisticsService();
        var probs = stats.ConvertChisToProbs(new[] { 10.0, 12.0, 15.0 });
        probs.Sum().Should().BeApproximately(1.0, 1e-9);
    }
}
