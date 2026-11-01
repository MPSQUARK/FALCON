using Falcon.Core.Numerics;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class SimdSpanOpsTests
{
    [Fact]
    public void AddWeighted_MatchesScalarReference()
    {
        var dest = new float[] { 1, 2, 3, 4, 5, 6, 7, 8 };
        var src = new float[] { 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f };
        var expected = (float[])dest.Clone();
        for (int i = 0; i < expected.Length; i++)
            expected[i] += 2f * src[i];

        SimdSpanOps.AddWeighted(dest, src, 2f);
        dest.Should().BeEquivalentTo(expected);
    }

    [Fact]
    public void ComputeChiRow_MatchesScalarReference()
    {
        var data = new float[] { 1, 2, 3, 4, 5, 6, 7, 8 };
        var error = new float[] { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f };
        var model = new float[] { 1.1f, 1.9f, 3.2f, 3.8f, 5.5f, 5.7f, 7.1f, 8.4f };
        var dest = new float[model.Length];
        var expected = new float[model.Length];
        for (int i = 0; i < model.Length; i++)
            expected[i] = (model[i] - data[i]) / error[i];

        SimdSpanOps.ComputeChiRow(dest, model, data, error);
        for (int i = 0; i < expected.Length; i++)
            dest[i].Should().BeApproximately(expected[i], 1e-5f);
    }

    [Fact]
    public void SampleMeanStd_MatchesScalarReference()
    {
        var values = new float[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        var (mean, std) = SimdSpanOps.SampleMeanStd(values);
        mean.Should().BeApproximately(5.5, 1e-6);
        std.Should().BeApproximately(3.02765, 1e-4);
    }

    [Fact]
    public void BlendWeightedAccumulate_MatchesScalarReference()
    {
        var dest = new double[] { 1, 2, 3, 4 };
        var src = new double[] { 10, 20, 30, 40 };
        var expected = new double[] { 3, 6, 9, 12 };
        SimdSpanOps.BlendWeightedAccumulate(dest, 0.2, src);
        dest.Should().BeEquivalentTo(expected);
    }

    [Fact]
    public void Sum_LargeN_KahanBeatsNaiveFloatAccumulation()
    {
        const int n = 100_000;
        var values = new float[n];
        for (int i = 0; i < n; i++)
            values[i] = 0.1f;

        float naive = 0f;
        for (int i = 0; i < n; i++)
            naive += values[i];

        double accurate = SimdSpanOps.Sum(values);
        accurate.Should().BeApproximately(n * 0.1, 1e-3);
        Math.Abs(naive - n * 0.1).Should().BeGreaterThan(1e-3, "sanity: naive f32 sum should drift for large N");
    }

    [Fact]
    public void SumSquares_LargeN_IsStable()
    {
        const int n = 8_192;
        var values = new float[n];
        for (int i = 0; i < n; i++)
            values[i] = 0.25f;

        SimdSpanOps.SumSquares(values).Should().BeApproximately(n * 0.25 * 0.25, n * 1e-4);
    }

    [Fact]
    public void AddWeighted_ManyTerms_UsesDoubleIntermediate()
    {
        const int n = 256;
        var dest = new float[n];
        var src = new float[n];
        Array.Fill(src, 1f);
        for (int m = 0; m < 200; m++)
            SimdSpanOps.AddWeighted(dest, src, 0.01f);

        for (int i = 0; i < n; i++)
            dest[i].Should().BeApproximately(2f, 1e-4f);
    }
}
