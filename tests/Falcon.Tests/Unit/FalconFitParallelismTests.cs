using Falcon.Core.Pipeline;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class FalconFitParallelismTests
{
    [Theory]
    [InlineData(52, 0, -1)]
    [InlineData(10, 4, 4)]
    [InlineData(3, 16, 3)]
    public void ResolveBatchParallelism_MatchesLegacyParallelFor(int files, int configured, int expected)
    {
        FalconFitParallelism.ResolveBatchParallelism(files, configured).Should().Be(expected);
    }

    [Fact]
    public void BatchParallelOptions_Default_IsUnlimitedLikeMasterProgram()
    {
        var options = FalconFitParallelism.BatchParallelOptions(52, configuredMax: 0);
        options.MaxDegreeOfParallelism.Should().Be(-1);
    }

    [Fact]
    public void DescribeParallelism_Unlimited_MentionsSpectrumCount()
    {
        FalconFitParallelism.DescribeParallelism(52, 0).Should().Contain("52");
    }
}
