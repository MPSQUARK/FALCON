using BAVCL;
using Falcon.Core.Abstractions;
using Falcon.Core.Diagnostics;
using Falcon.Core.Numerics;
using Falcon.Core.Services.Statistics;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;

namespace Falcon.Core.Services.Dust;

public sealed class DustService(GPU gpu) : IDustService
{
    readonly GPU _gpu = gpu;

    public double[] Hpf(ReadOnlySpan<double> flux, int windowSize = 20, int wStart = 20) =>
        GpuFftService.Hpf(_gpu, flux, windowSize, wStart);

    public double[] Unred(ReadOnlySpan<double> wavelengthAngstrom, double ebv, string law = "calzetti")
    {
        var wave = wavelengthAngstrom.ToArray();
        return law.ToLowerInvariant() switch
        {
            "allen" => DustAllen(ebv, wave),
            "prevot" => DustPrevot(ebv, wave),
            _ => DustCalzetti(ebv, wave)
        };
    }

    /// <inheritdoc />
    public double[] FitzpatrickUnred(ReadOnlySpan<double> wavelengthAngstrom, double ebv, double rV = 3.1)
    {
        var wave = wavelengthAngstrom.ToArray();
        int n = wave.Length;
        var x = new double[n];
        for (int i = 0; i < n; i++)
            x[i] = 10000.0 / wave[i];

        double x0 = 4.596;
        double gamma = 0.99;
        double c3 = 3.23;
        double c4 = 0.41;
        double c2 = -0.824 + 4.717 / rV;
        double c1 = 2.030 - 3.007 * c2;

        double xcutuv = 10000.0 / 2700.0;
        double[] xspluv = [10000.0 / 2700.0, 10000.0 / 2600.0];

        var curve = new double[n];
        var xuvList = new List<double>(xspluv);
        var yuvList = new List<double>();

        foreach (double xv in xspluv)
            yuvList.Add(FitzpatrickY(xv, x0, gamma, c3, c4, c1, c2, rV));

        for (int i = 0; i < n; i++)
        {
            if (x[i] >= xcutuv)
            {
                xuvList.Add(x[i]);
                yuvList.Add(FitzpatrickY(x[i], x0, gamma, c3, c4, c1, c2, rV));
                curve[i] = yuvList[^1];
            }
        }

        double[] xsplopir =
        [
            0,
            10000.0 / 26500.0, 10000.0 / 12200.0, 10000.0 / 6000.0,
            10000.0 / 5470.0, 10000.0 / 4670.0, 10000.0 / 4110.0,
            xspluv[0], xspluv[1]
        ];
        double[] ysplir = [0.0, 0.26469, 0.82925];
        for (int i = 0; i < ysplir.Length; i++)
            ysplir[i] *= rV / 3.1;
        double[] ysplop =
        [
            PolyVal([2.13572e-04, 1.00270, -4.22809e-01], rV),
            PolyVal([-7.35778e-05, 1.00216, -5.13540e-02], rV),
            PolyVal([-3.32598e-05, 1.00184, 7.00127e-01], rV),
            PolyVal([-4.45636e-05, 7.97809e-04, -5.46959e-03, 1.01707, 1.19456], rV)
        ];
        var ysplopir = ysplir.Concat(ysplop).Concat(yuvList.Take(2)).ToArray();
        var spline = CubicSpline.InterpolateNaturalSorted(xsplopir, ysplopir);

        for (int i = 0; i < n; i++)
        {
            if (x[i] < xcutuv)
                curve[i] = spline.Interpolate(x[i]);
        }

        var atten = new double[n];
        for (int i = 0; i < n; i++)
            atten[i] = Math.Pow(10, 0.4 * curve[i] * ebv);
        return atten;
    }

    static double FitzpatrickY(double xv, double x0, double gamma, double c3, double c4, double c1, double c2, double rV)
    {
        double y = c1 + c2 * xv;
        y += c3 * xv * xv / (Math.Pow(xv * xv - x0 * x0, 2) + Math.Pow(xv * gamma, 2));
        double bump = Math.Max(xv, 5.9) - 5.9;
        y += c4 * (0.5392 * bump * bump + 0.05644 * bump * bump * bump);
        return y + rV;
    }

    static double PolyVal(double[] coeffsHighToLow, double x)
    {
        double y = 0;
        foreach (double c in coeffsHighToLow)
            y = y * x + c;
        return y;
    }

    public double GetSfdEbv(double longitudeDeg, double latitudeDeg, float[,] ngp, float[,] sgp)
    {
        double l = longitudeDeg * Math.PI / 180;
        double b = latitudeDeg * Math.PI / 180;
        var map = b >= 0 ? ngp : sgp;
        int npix = map.GetLength(0);
        double theta = (Math.PI / 2) - b;
        double phi = l;
        if (phi < 0) phi += 2 * Math.PI;
        double x = phi / (2 * Math.PI) * npix;
        double y = theta / Math.PI * npix;
        int ix = (int)Math.Floor(x) % npix;
        int iy = Math.Clamp((int)Math.Floor(y), 0, npix - 1);
        return map[iy, ix];
    }

    public DustAttenuationResult DetermineAttenuation(DustAttenuationInput input, ISpectralFitter fitter)
    {
        using var _ = FitTrace.BeginSpan("dust_fit");
        var hpfData = Hpf(input.DataFlux);
        int nModels = input.ModelFlux.Length;
        var hpfModels = HpfModelsBatch(input.ModelFlux);

        for (int i = 0; i < hpfData.Length; i++)
            if (double.IsNaN(hpfData[i]) || double.IsInfinity(hpfData[i])) hpfData[i] = 0;

        var hpfError = new double[input.ErrorFlux.Length];
        double medErr = Median(input.ErrorFlux);
        double medData = Median(input.DataFlux);
        double medHpf = Median(hpfData);
        for (int i = 0; i < hpfError.Length; i++)
            hpfError[i] = medErr / medData * medHpf;

        var stats = new StatisticsService();
        var (normHpfModels, _) = stats.NormaliseSpec(hpfData, hpfModels);

        var flatModels = new double[nModels * input.Wavelength.Length];
        for (int m = 0; m < nModels; m++)
            Array.Copy(normHpfModels[m], 0, flatModels, m * input.Wavelength.Length, input.Wavelength.Length);

        var fitInput = new SpectralFitInput
        {
            Wavelength = input.Wavelength,
            Data = hpfData,
            Error = hpfError,
            Models = flatModels
        };
        var fit = fitter.Fit(fitInput);
        int best = fit.Solutions.OrderBy(s => s.ChiSquared).First().Index;
        var bestSol = fit.Solutions.First(s => s.Index == best);
        var bestHpf = Blend(bestSol.Weights, normHpfModels);
        var bestFull = Blend(bestSol.Weights, input.ModelFlux);

        var fine = new double[input.Wavelength.Length];
        for (int i = 0; i < fine.Length; i++)
        {
            fine[i] = (input.DataFlux[i] / (bestFull[i] == 0 ? 1e-30 : bestFull[i]))
                - (hpfData[i] / (bestHpf[i] == 0 ? 1e-30 : bestHpf[i])) + 1;
            if (double.IsNaN(fine[i]) || double.IsInfinity(fine[i])) fine[i] = 1;
        }
        fine = CurveSmoother(input.Wavelength, fine, input.SmoothingLength);
        double medFine = Median(fine);
        for (int i = 0; i < fine.Length; i++) fine[i] /= medFine;

        double bestEbv = 0;
        double bestChi = double.MaxValue;
        for (int ei = 0; ei < input.NumDustVals; ei++)
        {
            double e = ei / (input.MaxEbv * input.NumDustVals);
            var laws = Unred(input.Wavelength, e, input.DustLaw);
            double medLaw = Median(laws);
            for (int i = 0; i < laws.Length; i++) laws[i] /= medLaw;
            double chi = 0;
            for (int i = 0; i < fine.Length; i++)
            {
                double d = fine[i] - laws[i];
                chi += d * d;
            }
            if (chi < bestChi) { bestChi = chi; bestEbv = e; }
        }

        var attenuation = Unred(input.Wavelength, bestEbv, input.DustLaw);
        double avAge = Dot(bestSol.Weights, input.Age);
        double avMetal = Dot(bestSol.Weights, input.Metal);

        return new DustAttenuationResult
        {
            Ebv = bestEbv,
            AttenuationCurve = attenuation,
            AvAgeHpf = avAge,
            AvMetalHpf = avMetal
        };
    }

    private static double[] Blend(double[] weights, double[][] models)
    {
        int len = models[0].Length;
        var result = new double[len];
        for (int m = 0; m < weights.Length; m++)
            SimdSpanOps.BlendWeightedAccumulate(result, weights[m], models[m]);
        return result;
    }

    private static double Dot(double[] w, double[] v) => SimdSpanOps.Dot(w, v);

    double[][] HpfModelsBatch(double[][] models)
    {
        var hpfModels = new double[models.Length][];
        for (int m = 0; m < models.Length; m++)
            hpfModels[m] = Hpf(models[m]);
        return hpfModels;
    }

    private static double[] CurveSmoother(double[] x, double[] y, float smoothingLength)
    {
        var yOut = new double[y.Length];
        var scratch = new double[y.Length];
        int lo = 0;
        int hi = 0;
        for (int w = 0; w < x.Length; w++)
        {
            double left = x[w] - smoothingLength;
            double right = x[w] + smoothingLength;
            while (lo < x.Length && x[lo] < left) lo++;
            hi = Math.Max(hi, lo);
            while (hi < x.Length && x[hi] <= right) hi++;
            int count = hi - lo;
            if (count <= 0)
            {
                yOut[w] = y[w];
                continue;
            }
            Array.Copy(y, lo, scratch, 0, count);
            Array.Sort(scratch, 0, count);
            yOut[w] = count % 2 == 0
                ? (scratch[count / 2 - 1] + scratch[count / 2]) / 2
                : scratch[count / 2];
        }
        return yOut;
    }

    private static double[] DustCalzetti(double ebv, double[] lam)
    {
        var k = new double[lam.Length];
        for (int i = 0; i < lam.Length; i++)
        {
            double l = lam[i] / 10000.0;
            if (l >= 0.63 && l <= 2.2) k[i] = 2.659 * (-1.857 + 1.040 / l) + 4.05;
            else if (l < 0.63) k[i] = 2.659 * (-2.156 + 1.509 / l - 0.198 / (l * l) + 0.011 / (l * l * l)) + 4.05;
            else k[i] = 0;
        }
        return k.Select(v => Math.Pow(10, -0.4 * ebv * v)).ToArray();
    }

    private static double[] DustAllen(double ebv, double[] lam)
    {
        double[] wave = [1000, 1110, 1250, 1430, 1670, 2000, 2220, 2500, 2860, 3330, 3650, 4000, 4400, 5000, 5530, 6700, 9000, 10000];
        double[] allenK = [4.20, 3.70, 3.30, 3.00, 2.70, 2.80, 2.90, 2.30, 1.97, 1.69, 1.58, 1.45, 1.32, 1.13, 1.00, 0.74, 0.46, 0.38];
        for (int i = 0; i < allenK.Length; i++) allenK[i] *= 3.1;
        return lam.Select(l => Math.Pow(10, -0.4 * ebv * InterpolateTable(wave, allenK, l))).ToArray();
    }

    private static double[] DustPrevot(double ebv, double[] lam)
    {
        double[] wave = [1275, 1330, 1385, 1435, 1490, 1545, 1595, 1647, 1700, 1755, 1810, 1860, 1910, 2000, 2115, 2220, 2335, 2445, 2550, 2665, 2778, 2890, 2995, 3105, 3704, 4255, 5291, 10000];
        double[] prevotK = [13.54, 12.52, 11.51, 10.80, 9.84, 9.28, 9.06, 8.49, 8.01, 7.71, 7.17, 6.90, 6.76, 6.38, 5.85, 5.30, 4.53, 4.24, 3.91, 3.49, 3.15, 3.00, 2.65, 2.29, 1.81, 1.00, 0.74, 0.00];
        for (int i = 0; i < prevotK.Length; i++) prevotK[i] *= 2.72;
        return lam.Select(l => Math.Pow(10, -0.4 * ebv * InterpolateTable(wave, prevotK, l))).ToArray();
    }

    private static double InterpolateTable(double[] wave, double[] k, double lam)
    {
        if (lam <= wave[0]) return k[0];
        if (lam >= wave[^1]) return k[^1];
        int idx = Array.BinarySearch(wave, lam);
        if (idx >= 0) return k[idx];
        idx = ~idx;
        double t = (lam - wave[idx - 1]) / (wave[idx] - wave[idx - 1]);
        return k[idx - 1] * (1 - t) + k[idx] * t;
    }

    private static double Median(ReadOnlySpan<double> data)
    {
        var sorted = data.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2 : sorted[n / 2];
    }

    private static double Median(IEnumerable<double> data) => Median(data.ToArray());

}
