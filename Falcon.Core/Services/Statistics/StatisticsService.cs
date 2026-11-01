using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.Statistics;

public sealed class StatisticsService : IStatisticsService
{
    public (double[][] NormalizedModels, double[] MassFactors) NormaliseSpec(ReadOnlySpan<double> dataFlux, double[][] modelFlux)
    {
        double dataNorm = Median(dataFlux);
        int numMods = modelFlux.Length;
        var normed = new double[numMods][];
        var massFactors = new double[numMods];
        for (int m = 0; m < numMods; m++)
        {
            double modelNorm = Median(modelFlux[m]);
            massFactors[m] = modelNorm == 0 ? 1 : dataNorm / modelNorm;
            normed[m] = modelFlux[m].Select(v => modelNorm == 0 ? v : v / modelNorm * dataNorm).ToArray();
        }
        return (normed, massFactors);
    }

    public (double[] MatchedWave, double[] MatchedData, double[] MatchedError, double[][] MatchedModels) MatchDataModels(
        double[] dataWave, double[] dataFlux, bool[] dataFlags, double[] errorFlux,
        double[] modelWave, double[][] modelFlux, double minWave, double maxWave)
    {
        var flaggedWave = dataWave.Where((_, i) => dataFlags[i]).ToArray();
        if (flaggedWave.Length == 0) flaggedWave = dataWave;

        double minW = Math.Max(flaggedWave.Min(), Math.Max(modelWave.Min(), minWave));
        double maxW = Math.Min(flaggedWave.Max(), Math.Min(modelWave.Max(), maxWave));

        int numModels = modelFlux.Length;
        var locModel = modelWave.Select(w => w >= minW && w <= maxW).ToArray();
        int numMod = locModel.Count(x => x);
        if (numMod == 0) throw new ArgumentException("Wavelength range outside model coverage.");
        var trimmedModelWave = modelWave.Where((w, i) => locModel[i]).ToArray();
        var trimmedModelFlux = new double[numModels][];
        for (int m = 0; m < numModels; m++)
            trimmedModelFlux[m] = modelFlux[m].Where((_, i) => locModel[i]).ToArray();

        var locData = dataWave.Select(w => w >= minW && w <= maxW).ToArray();
        int numDat = locData.Count(x => x);
        if (numDat == 0) throw new ArgumentException("Wavelength range outside data coverage.");
        var trimmedDataWave = dataWave.Where((_, i) => locData[i]).ToArray();
        var trimmedDataFlux = dataFlux.Where((_, i) => locData[i]).ToArray();
        var trimmedError = errorFlux.Where((_, i) => locData[i]).ToArray();

        if (numMod >= numDat)
        {
            var bisect = StatisticsMath.BisectArray(trimmedDataWave);
            for (int i = 0; i < bisect.Length; i++)
                bisect[i] += trimmedDataWave.Min() * 1e-10;
            var matchedModel = new double[numModels][];
            for (int m = 0; m < numModels; m++)
                matchedModel[m] = StatisticsMath.IntegrateToBins(trimmedModelWave, trimmedModelFlux[m], bisect);
            return (trimmedDataWave[1..^1], trimmedDataFlux[1..^1], trimmedError[1..^1], matchedModel);
        }

        var bisectModel = StatisticsMath.BisectArray(trimmedModelWave);
        for (int i = 0; i < bisectModel.Length; i++)
            bisectModel[i] += trimmedModelWave.Min() * 1e-10;
        var matchedData = StatisticsMath.IntegrateToBins(trimmedDataWave, trimmedDataFlux, bisectModel);
        var matchedError = StatisticsMath.IntegrateToBins(trimmedDataWave, trimmedError, bisectModel);
        var matchedModels = new double[numModels][];
        for (int m = 0; m < numModels; m++)
            matchedModels[m] = trimmedModelFlux[m][1..^1];
        return (trimmedModelWave[1..^1], matchedData, matchedError, matchedModels);
    }

    public double[] ConvertChisToProbs(ReadOnlySpan<double> chiSquared, int dof) =>
        StatisticsMath.ConvertChisToProbs(chiSquared, dof);

    public Dictionary<string, double> CalculateAveragesPdf(
        double[] probs, double[][] lightWeights, double[][] massWeights, double[][] unnormMass,
        double[] age, double[] metal, double fluxUnits, double distLumCm, int sampling) =>
        StatisticsMath.CalculateAveragesPdf(probs, lightWeights, massWeights, unnormMass, age, metal, sampling, distLumCm, fluxUnits);

    public (double[][] UnnormMass, double[][] MassWeights) LightWeightsToMass(double[][] lightWeights, double[] massFactors) =>
        StatisticsMath.LightWeightsToMass(lightWeights, massFactors);

    public double[] LightWeightsToMassRow(double[] lightWeights, double[] massFactors) =>
        StatisticsMath.LightWeightsToMassRow(lightWeights, massFactors);

    private static double Median(ReadOnlySpan<double> data)
    {
        var sorted = data.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2 : sorted[n / 2];
    }
}
