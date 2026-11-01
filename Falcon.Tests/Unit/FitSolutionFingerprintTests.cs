using Falcon.Core.Services.Fitting;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class FitSolutionFingerprintTests
{
    [Fact]
    public void FromWeights_TreatsNearPureSingleModelAsEquivalent()
    {
        var pure = new double[200];
        pure[188] = 1.0;

        var nearPure = new double[200];
        nearPure[188] = 0.99995;
        nearPure[187] = 0.00005;

        FitSolutionFingerprint.FromWeights(pure)
            .Should()
            .Be(FitSolutionFingerprint.FromWeights(nearPure));
    }

    [Fact]
    public void FromWeights_DistinguishesDifferentDominantModels()
    {
        var a = new double[200];
        a[188] = 1.0;
        var b = new double[200];
        b[189] = 1.0;

        FitSolutionFingerprint.FromWeights(a)
            .Should()
            .NotBe(FitSolutionFingerprint.FromWeights(b));
    }

    [Fact]
    public void FromWeights_DistinguishesMeaningfulBlends()
    {
        var a = new double[4];
        a[0] = 0.7;
        a[1] = 0.3;

        var b = new double[4];
        b[0] = 0.6;
        b[1] = 0.4;

        FitSolutionFingerprint.FromWeights(a)
            .Should()
            .NotBe(FitSolutionFingerprint.FromWeights(b));
    }
}
