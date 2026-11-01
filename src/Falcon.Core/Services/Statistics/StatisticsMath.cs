using MathNet.Numerics.Distributions;

namespace Falcon.Core.Services.Statistics;

/// <summary>Static firefly_library equivalents exposed for unit testing.</summary>
public static class StatisticsMath
{
    public static double[] BisectArray(ReadOnlySpan<double> array)
    {
        var bins = new List<double> { array[0] };
        for (int i = 1; i < array.Length; i++)
            bins.Add((array[i] + array[i - 1]) / 2);
        bins.Add(array[^1]);
        return bins.ToArray();
    }

    public static (double[] PropertyPdf, double[] ProbPdf) MaxPdf(ReadOnlySpan<double> probs, ReadOnlySpan<double> property, int sampling)
    {
        double lower = property.ToArray().Min();
        double upper = property.ToArray().Max();
        double errorInterval = Math.Round(upper, 2) - Math.Round(lower, 2);

        if (Math.Round(upper, 2) == Math.Round(lower, 2) || errorInterval <= Math.Abs(upper / 100.0 * 3))
            return (property.ToArray(), Enumerable.Repeat(1.0 / probs.Length, probs.Length).ToArray());

        var propertyPdfInt = new List<double>();
        double step = (upper - lower) / sampling;
        for (double v = lower; v <= upper * 1.001; v += step)
            propertyPdfInt.Add(v + (upper - lower) * 0.000001);

        var probPdf = new double[propertyPdfInt.Count];
        for (int p = 0; p < propertyPdfInt.Count - 1; p++)
        {
            double lo = propertyPdfInt[p];
            double hi = propertyPdfInt[p + 1];
            var matches = new List<double>();
            for (int i = 0; i < property.Length; i++)
                if (property[i] <= hi && property[i] > lo)
                    matches.Add(probs[i]);
            probPdf[p] = matches.Count > 0 ? matches.Max() : 0;
        }

        var propertyPdf = BisectArray(propertyPdfInt.ToArray());
        double sum = probPdf.Take(probPdf.Length - 1).Sum();
        var normalized = probPdf.Take(probPdf.Length - 1).Select(v => v / (sum == 0 ? 1 : sum)).ToArray();
        return (propertyPdf, normalized);
    }

    public static double[] ConvertChisToProbs(ReadOnlySpan<double> chiSquared, int dof)
    {
        if (chiSquared.IsEmpty) return Array.Empty<double>();
        double minChi = chiSquared.ToArray().Min();
        if (minChi <= 0) minChi = 1e-30;
        var chi2 = new ChiSquared(dof);
        var scaled = chiSquared.ToArray().Select(c => c / minChi * dof).ToArray();
        var probs = scaled.Select(c => 1.0 - chi2.CumulativeDistribution(c)).ToArray();
        double sum = probs.Sum();
        return probs.Select(p => p / (sum == 0 ? 1 : sum)).ToArray();
    }

    public static (double[][] UnnormMass, double[][] MassWeights) LightWeightsToMass(double[][] lightWeights, double[] massFactors)
    {
        int n = lightWeights.Length;
        var unnorm = new double[n][];
        var mass = new double[n][];
        for (int w = 0; w < n; w++)
        {
            int cols = lightWeights[w].Length;
            unnorm[w] = new double[cols];
            mass[w] = new double[cols];
            for (int c = 0; c < cols; c++)
                unnorm[w][c] = lightWeights[w][c] * massFactors[c];
            double rowSum = unnorm[w].Sum();
            for (int c = 0; c < cols; c++)
                mass[w][c] = rowSum == 0 ? 0 : unnorm[w][c] / rowSum;
        }
        return (unnorm, mass);
    }

    public static double[] LightWeightsToMassRow(double[] lightWeights, double[] massFactors)
    {
        var unnorm = lightWeights.Zip(massFactors, (w, f) => w * f).ToArray();
        double sum = unnorm.Sum();
        return sum == 0 ? unnorm : unnorm.Select(u => u / sum).ToArray();
    }

    public static (double UnnormMass, double[] MassWeights) LightWeightsToMassSingle(double[] lightWeights, double[] massFactors)
    {
        var unnorm = lightWeights.Zip(massFactors, (w, f) => w * f).ToArray();
        double sum = unnorm.Sum();
        var mass = sum == 0 ? unnorm : unnorm.Select(u => u / sum).ToArray();
        return (sum, mass);
    }

    public static int FindClosest(ReadOnlySpan<double> sortedA, double target)
    {
        int idx = Array.BinarySearch(sortedA.ToArray(), target);
        if (idx < 0) idx = ~idx;
        idx = Math.Clamp(idx, 1, sortedA.Length - 1);
        double left = sortedA[idx - 1];
        double right = sortedA[idx];
        idx -= target - left < right - target ? 1 : 0;
        return idx;
    }

    public static double[] AveragesAndErrors(ReadOnlySpan<double> probs, ReadOnlySpan<double> prop, int sampling)
    {
        var valid = prop.ToArray().Where(v => !double.IsNaN(v)).ToArray();
        if (probs.Length <= 1 || valid.Length <= 1)
            return [999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0];

        var (xdf, y) = MaxPdf(probs, prop, sampling);
        var cdf = new double[y.Length];
        for (int m = 0; m < y.Length; m++)
            cdf[m] = y.Take(m).Sum();
        double maxCdf = cdf.Max();
        if (maxCdf > 0)
            for (int m = 0; m < cdf.Length; m++) cdf[m] /= maxCdf;

        double dx = xdf.Length > 1 ? xdf[1] - xdf[0] : 1;
        var area = y.Select(v => v * dx).ToArray();
        double areaSum = area.Sum();
        var areaNorm = area.Select(a => areaSum == 0 ? 0 : a / areaSum).ToArray();
        var indx = areaNorm.Select((v, i) => (v, i)).OrderByDescending(t => t.v).Select(t => t.i).ToArray();
        var desc = indx.Select(i => areaNorm[i]).ToArray();
        var cdfProb = new double[desc.Length];
        for (int m = 0; m < desc.Length; m++)
            cdfProb[m] = desc.Take(m).Sum();

        double[] avSigs = [0.6827, 0.9545, 0.9973];
        var indexClose = avSigs.Select(s => FindClosest(cdfProb, s)).ToArray();

        double best = xdf[indx[0]];
        double u1 = indx.Take(indexClose[0] + 1).Max(i => xdf[i]);
        double l1 = indx.Take(indexClose[0] + 1).Min(i => xdf[i]);
        double u2 = indx.Take(indexClose[1] + 1).Max(i => xdf[i]);
        double l2 = indx.Take(indexClose[1] + 1).Min(i => xdf[i]);
        double u3 = indx.Take(indexClose[2] + 1).Max(i => xdf[i]);
        double l3 = indx.Take(indexClose[2] + 1).Min(i => xdf[i]);
        return [best, u1, l1, u2, l2, u3, l3];
    }

    public static Dictionary<string, double> CalculateAveragesPdf(
        double[] probs, double[][] lightWeights, double[][] massWeights, double[][] unnormMass,
        double[] age, double[] metal, int sampling, double distLumCm, double fluxUnits)
    {
        var av = new Dictionary<string, double>();
        int nSolutions = probs.Length;
        int nModels = age.Length;

        var lightAge = new double[nSolutions];
        var lightMetal = new double[nSolutions];
        var massAge = new double[nSolutions];
        var massMetal = new double[nSolutions];
        var totMass = new double[nSolutions];
        double conversion = fluxUnits * 4 * Math.PI * distLumCm * distLumCm;

        for (int s = 0; s < nSolutions; s++)
        {
            lightAge[s] = Dot(lightWeights[s], age);
            lightMetal[s] = Dot(lightWeights[s], metal);
            massAge[s] = Dot(massWeights[s], age);
            massMetal[s] = Dot(massWeights[s], metal);
            totMass[s] = unnormMass[s].Sum() * conversion;
        }

        SetAverages(av, "light_age", probs, lightAge, sampling);
        SetAverages(av, "light_metal", probs, lightMetal, sampling);
        SetAverages(av, "mass_age", probs, massAge, sampling);
        SetAverages(av, "mass_metal", probs, massMetal, sampling);
        SetAverages(av, "stellar_mass", probs, totMass, sampling);
        return av;
    }

    private static void SetAverages(Dictionary<string, double> av, string key, double[] probs, double[] values, int sampling)
    {
        var stats = AveragesAndErrors(probs, values, sampling);
        av[key] = stats[0];
        av[$"{key}_1_sig_plus"] = stats[1];
        av[$"{key}_1_sig_minus"] = stats[2];
        av[$"{key}_2_sig_plus"] = stats[3];
        av[$"{key}_2_sig_minus"] = stats[4];
        av[$"{key}_3_sig_plus"] = stats[5];
        av[$"{key}_3_sig_minus"] = stats[6];
    }

    private static double Dot(double[] weights, double[] values)
    {
        double s = 0;
        for (int i = 0; i < weights.Length; i++) s += weights[i] * values[i];
        return s;
    }

    public static double[] IntegrateToBins(double[] wave, double[] flux, double[] bisect)
    {
        var combinedWave = wave.Concat(bisect).OrderBy(x => x).ToArray();
        var interp = combinedWave.Select(x => Interpolate(wave, flux, x)).ToArray();

        var output = new double[bisect.Length - 1];
        for (int l = 0; l < output.Length; l++)
        {
            int start = Array.BinarySearch(combinedWave, bisect[l]);
            if (start < 0) start = ~start;
            int end = Array.BinarySearch(combinedWave, bisect[l + 1]);
            if (end < 0) end = ~end;
            if (end >= combinedWave.Length) { output[l] = l > 0 ? output[l - 1] : 0; continue; }
            double integral = 0;
            for (int i = start; i < end; i++)
                integral += interp[i] * (combinedWave[i + 1] - combinedWave[i]);
            output[l] = integral / (bisect[l + 1] - bisect[l]);
        }
        return output;
    }

    public static double Interpolate(double[] wave, double[] flux, double x)
    {
        int idx = Array.BinarySearch(wave, x);
        if (idx >= 0) return flux[idx];
        idx = ~idx;
        if (idx <= 0) return flux[0];
        if (idx >= wave.Length) return flux[^1];
        double t = (x - wave[idx - 1]) / (wave[idx] - wave[idx - 1]);
        return flux[idx - 1] * (1 - t) + flux[idx] * t;
    }
}
