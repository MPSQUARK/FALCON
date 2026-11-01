using Falcon.Core.Services.Statistics;
using Falcon.Tests.Helpers;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class MatchDataModelsTests
{
    private readonly StatisticsService _stats = new();

    [Fact]
    public void MatchDataModels_MatchesGoldenFixture()
    {
        using var golden = FireflyReferenceLoader.LoadGolden("golden_match_data_models.json");
        var dataWave = new[] { 4000.0, 4010.0, 4020.0, 4030.0, 4040.0 };
        var dataFlux = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
        var err = new[] { 0.1, 0.1, 0.1, 0.1, 0.1 };
        var flags = new[] { true, true, true, true, true };
        var modelWave = new[] { 4005.0, 4015.0, 4025.0 };
        var modelFlux = new[] { new[] { 1.0, 2.0, 3.0 }, new[] { 2.0, 3.0, 4.0 } };

        var match = _stats.MatchDataModels(dataWave, dataFlux, flags, err, modelWave, modelFlux, 4000, 4040);
        GoldenTestHelpers.AssertArraysEqual(match.MatchedWave, FireflyReferenceLoader.GetDoubleArray(golden, "matched_wave"));
        GoldenTestHelpers.AssertArraysEqual(match.MatchedData, FireflyReferenceLoader.GetDoubleArray(golden, "matched_data"));
        GoldenTestHelpers.AssertArraysEqual(match.MatchedModels[0], FireflyReferenceLoader.GetDoubleArray(golden, "matched_model_0"));
    }

    [Fact]
    public void MatchDataModels_DataDenserThanModels_RegridsOntoModelAxis()
    {
        var dataWave = Enumerable.Range(0, 10).Select(i => 4000.0 + i * 10).ToArray();
        var dataFlux = dataWave.Select(w => w * 1e-4).ToArray();
        var err = Enumerable.Repeat(0.01, 10).ToArray();
        var flags = Enumerable.Repeat(true, 10).ToArray();
        var modelWave = new[] { 4005.0, 4025.0, 4045.0, 4065.0, 4085.0 };
        var modelFlux = new[] { modelWave.Select(w => w * 1e-4).ToArray() };

        var match = _stats.MatchDataModels(dataWave, dataFlux, flags, err, modelWave, modelFlux, 4000, 4100);
        match.MatchedWave.Length.Should().Be(modelWave.Length - 2);
        match.MatchedModels[0].Length.Should().Be(match.MatchedWave.Length);
    }

    [Fact]
    public void MatchDataModels_ModelDenserThanData_RegridsOntoDataAxis()
    {
        var modelWave = Enumerable.Range(0, 10).Select(i => 4000.0 + i * 10).ToArray();
        var modelFlux = new[] { modelWave.Select(w => w * 1e-4).ToArray() };
        var dataWave = new[] { 4005.0, 4025.0, 4045.0, 4065.0, 4085.0 };
        var dataFlux = dataWave.Select(w => w * 1e-4).ToArray();
        var err = Enumerable.Repeat(0.01, 5).ToArray();
        var flags = Enumerable.Repeat(true, 5).ToArray();

        var match = _stats.MatchDataModels(dataWave, dataFlux, flags, err, modelWave, modelFlux, 4000, 4100);
        match.MatchedWave.Length.Should().Be(dataWave.Length - 2);
        match.MatchedData.Length.Should().Be(match.MatchedWave.Length);
    }

    [Fact]
    public void MatchDataModels_EmptyFlaggedRange_UsesFullDataWave()
    {
        var dataWave = new[] { 4000.0, 4010.0, 4020.0, 4030.0, 4040.0 };
        var dataFlux = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
        var err = Enumerable.Repeat(0.1, 5).ToArray();
        var flags = new[] { false, false, false, false, false };
        var modelWave = new[] { 4005.0, 4015.0, 4025.0 };
        var modelFlux = new[] { new[] { 1.0, 2.0, 3.0 } };

        var act = () => _stats.MatchDataModels(dataWave, dataFlux, flags, err, modelWave, modelFlux, 4000, 4040);
        act.Should().NotThrow();
    }

    [Fact]
    public void MatchDataModels_OutOfRangeWave_Throws()
    {
        var dataWave = new[] { 5000.0, 5010.0 };
        var dataFlux = new[] { 1.0, 2.0 };
        var err = new[] { 0.1, 0.1 };
        var flags = new[] { true, true };
        var modelWave = new[] { 4000.0, 4010.0 };
        var modelFlux = new[] { new[] { 1.0, 2.0 } };

        var act = () => _stats.MatchDataModels(dataWave, dataFlux, flags, err, modelWave, modelFlux, 4000, 4100);
        act.Should().Throw<ArgumentException>();
    }
}
