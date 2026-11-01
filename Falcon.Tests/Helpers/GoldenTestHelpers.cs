using FluentAssertions;

namespace Falcon.Tests.Helpers;

public static class GoldenTestHelpers
{
    public const double DefaultRtol = 1e-5;

    public static void AssertArraysEqual(IReadOnlyList<double> actual, IReadOnlyList<double> expected, double rtol = DefaultRtol)
    {
        actual.Count.Should().Be(expected.Count);
        for (int i = 0; i < actual.Count; i++)
            actual[i].Should().BeApproximately(expected[i], Math.Max(rtol, rtol * Math.Abs(expected[i]) + 1e-12));
    }

    public static void AssertSumToOne(double[] probs, double tolerance = 1e-10) =>
        probs.Sum().Should().BeApproximately(1.0, tolerance);
}
