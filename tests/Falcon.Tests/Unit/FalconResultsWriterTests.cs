using Falcon.Core.Abstractions;
using Falcon.Core.Services.IO;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class FalconResultsWriterTests
{
    [Fact]
    public void WriteAndLoad_RoundTripsTopSolutions()
    {
        var dir = Path.Combine(Path.GetTempPath(), "falcon-results-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        var resultsPath = Path.Combine(dir, "spec.results.json");

        var output = new FireflyFitOutput
        {
            SourceFile = "assets/spec.fits",
            Converged = true,
            BestEbv = 0.05,
            MinChiSquared = 10,
            BestLightWeights = [1],
            BestMassWeights = [1],
            Averages = new Dictionary<string, double> { ["light_age"] = 5e9 },
            MatchedWave = [4000, 4100, 4200],
            MatchedData = [1, 2, 1.5],
            MatchedError = [0.1, 0.1, 0.1],
            BestFitSpectrum = [1, 2, 1.5],
            OutputPath = Path.Combine(dir, "spec_firefly.fits"),
            ResultsPath = resultsPath,
            TopSolutions =
            [
                new FalconSolutionRecord
                {
                    Rank = 1,
                    ChiSquared = 10,
                    ChiSquaredPerDof = 3.33,
                    Probability = 0.8,
                    LightWeightedAgeGyr = 5,
                    LightWeightedMetal = 0.02,
                    ModelFlux = [1.01, 1.99, 1.48],
                    Components = [new FalconSspComponent { ModelIndex = 0, LightWeight = 1, AgeGyr = 5, Metal = -0.2 }]
                }
            ]
        };

        FalconResultsWriter.Write(output, redshift: 0.02f);
        File.Exists(resultsPath).Should().BeTrue();

        var doc = FalconResultsReader.Load(resultsPath);
        doc.TopSolutions.Should().HaveCount(1);
        doc.WavelengthAngstrom.Should().HaveCount(3);
        doc.DataFlux[1].Should().Be(2);
        doc.TopSolutions[0].ModelFlux[0].Should().BeApproximately(1.01, 1e-6);
    }

    [Fact]
    public void ComputeModelWaveRange_TrimsZeroModelEdges()
    {
        var wave = new[] { 3500.0, 3600.0, 3700.0, 3800.0, 3900.0 };
        var model = new[] { 0.0, 1.0, 2.0, 1.0, 0.0 };
        var (wMin, wMax) = FalconResultsWriter.ComputeModelWaveRange(wave, model);
        wMin.Should().Be(3600);
        wMax.Should().Be(3800);
    }
}
