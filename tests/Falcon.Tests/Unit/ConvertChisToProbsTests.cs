using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class ConvertChisToProbsTests
{
    private readonly StatisticsService _stats = new();

    [Fact]
    public void ConvertChisToProbs_SumsToOne()
    {
        var probs = _stats.ConvertChisToProbs([10.0, 12.0, 15.0, 20.0], dof: 4);
        GoldenTestHelpers.AssertSumToOne(probs);
    }

    [Fact]
    public void ConvertChisToProbs_EmptyInput_ReturnsEmpty()
    {
        StatisticsMath.ConvertChisToProbs([], dof: 3).Should().BeEmpty();
    }

    [Fact]
    public void ConvertChisToProbs_BestChiGetsHighestProbability()
    {
        var chi = new[] { 5.0, 10.0, 15.0 };
        var probs = _stats.ConvertChisToProbs(chi, dof: 2);
        probs[0].Should().BeGreaterThan(probs[1]);
        probs[1].Should().BeGreaterThan(probs[2]);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(5)]
    [InlineData(10)]
    public void ConvertChisToProbs_VariousDof_AllNormalize(int dof)
    {
        var chi = new[] { 1.0, 2.0, 3.0, 4.0 };
        var probs = _stats.ConvertChisToProbs(chi, dof);
        GoldenTestHelpers.AssertSumToOne(probs);
        probs.Should().OnlyContain(p => p >= 0);
    }

    [Theory]
    [InlineData(10.0, 10.0, 10.0)]
    [InlineData(1.0, 1.01, 1.02)]
    public void ConvertChisToProbs_ServiceMatchesStatisticsMath(double c0, double c1, double c2)
    {
        var chi = new[] { c0, c1, c2 };
        var fromService = _stats.ConvertChisToProbs(chi, 3);
        var fromMath = StatisticsMath.ConvertChisToProbs(chi, 3);
        GoldenTestHelpers.AssertArraysEqual(fromService, fromMath);
    }

    [Fact]
    public void ConvertChisToProbs_SingleValue_IsOne()
    {
        var probs = _stats.ConvertChisToProbs([42.0], dof: 5);
        probs.Should().ContainSingle().Which.Should().BeApproximately(1.0, 1e-9);
    }
}
