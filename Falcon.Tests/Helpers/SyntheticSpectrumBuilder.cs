namespace Falcon.Tests.Helpers;

public static class SyntheticSpectrumBuilder
{
    public static (double[] Wave, double[] Flux, double[] Error) LinearGrid(int n, double start = 4000, double step = 10)
    {
        var wave = Enumerable.Range(0, n).Select(i => start + i * step).Select(d => (double)d).ToArray();
        var flux = wave.Select(w => 1.0 + w * 1e-5).ToArray();
        var err = Enumerable.Repeat(0.01, n).Select(d => (double)d).ToArray();
        return (wave, flux, err);
    }

    public static double[][] ConstantModels(int nModels, int nWave, double value = 1.0)
    {
        var models = new double[nModels][];
        for (int m = 0; m < nModels; m++)
            models[m] = Enumerable.Repeat(value + m * 0.1, nWave).Select(d => (double)d).ToArray();
        return models;
    }
}
