using System.Numerics;
using Falcon.Core.Abstractions;
using Falcon.Core.Services.Statistics;

namespace Falcon.Core.Services.Dust;

public sealed class DustService : IDustService
{
    public double[] Hpf(ReadOnlySpan<double> flux, int windowSize = 20, int wStart = 20)
    {
        int d = flux.Length;
        if (wStart == 0 && windowSize == 0) wStart = Math.Max(1, d / 100);

        var h = FftForward(flux.ToArray());
        var window = new double[d];
        var unwindow = new double[d];
        window[0] = 1;
        unwindow[0] = 1;

        if (windowSize > 0)
        {
            double dw = windowSize;
            for (int i = 0; i < windowSize; i++)
            {
                window[wStart + i] = (i + 1.0) / dw;
                window[d - 1 - (wStart + windowSize - i)] = (dw - i) / dw;
            }
            for (int i = wStart + windowSize; i < d - (wStart + windowSize); i++)
                window[i] = 1;
        }
        else
        {
            for (int i = wStart; i < d - wStart; i++) window[i] = 1;
        }

        for (int i = 0; i < d; i++) unwindow[i] = 1 - window[i];
        unwindow[0] = 1;

        var hFiltered = new Complex[d];
        var unHFiltered = new Complex[d];
        for (int i = 0; i < d; i++)
        {
            hFiltered[i] = h[i] * window[i];
            unHFiltered[i] = h[i] * unwindow[i];
        }

        var res = FftInverseReal(hFiltered);
        var unres = FftInverseReal(unHFiltered);
        double medRes = Median(res);
        var output = new double[d];
        for (int i = 0; i < d; i++)
            output[i] = (1.0 + (res[i] - medRes) / (unres[i] == 0 ? 1e-30 : unres[i])) * medRes;
        return output;
    }

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
        var hpfData = Hpf(input.DataFlux);
        int nModels = input.ModelFlux.Length;
        var hpfModels = new double[nModels][];
        for (int m = 0; m < nModels; m++) hpfModels[m] = Hpf(input.ModelFlux[m]);

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
            for (int i = 0; i < len; i++)
                result[i] += weights[m] * models[m][i];
        return result;
    }

    private static double Dot(double[] w, double[] v)
    {
        double s = 0;
        for (int i = 0; i < w.Length; i++) s += w[i] * v[i];
        return s;
    }

    private static double[] CurveSmoother(double[] x, double[] y, float smoothingLength)
    {
        var yOut = new double[y.Length];
        for (int w = 0; w < x.Length; w++)
        {
            var vals = new List<double>();
            for (int i = 0; i < x.Length; i++)
                if (x[i] < x[w] + smoothingLength && x[i] > x[w] - smoothingLength)
                    vals.Add(y[i]);
            yOut[w] = vals.Count > 0 ? Median(vals) : y[w];
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
        // Simplified Allen 1976 MW curve via Calzetti fallback coefficients for parity scaffolding
        return DustCalzetti(ebv * 1.1, lam);
    }

    private static double[] DustPrevot(double ebv, double[] lam) => DustCalzetti(ebv * 0.9, lam);

    private static double Median(ReadOnlySpan<double> data)
    {
        var sorted = data.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2 : sorted[n / 2];
    }

    private static double Median(IEnumerable<double> data) => Median(data.ToArray());

    private static Complex[] FftForward(double[] input)
    {
        int n = input.Length;
        int pow2 = 1;
        while (pow2 < n) pow2 <<= 1;
        var padded = new Complex[pow2];
        for (int i = 0; i < n; i++) padded[i] = new Complex(input[i], 0);
        FftInPlace(padded, false);
        return padded;
    }

    private static double[] FftInverseReal(Complex[] spectrum)
    {
        var copy = (Complex[])spectrum.Clone();
        FftInPlace(copy, true);
        return copy.Select(c => c.Real / copy.Length).ToArray();
    }

    private static void FftInPlace(Complex[] buffer, bool inverse)
    {
        int n = buffer.Length;
        int j = 0;
        for (int i = 1; i < n; i++)
        {
            int bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1) j &= ~bit;
            j |= bit;
            if (i < j) (buffer[i], buffer[j]) = (buffer[j], buffer[i]);
        }
        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = 2 * Math.PI / len * (inverse ? -1 : 1);
            var wlen = new Complex(Math.Cos(ang), Math.Sin(ang));
            for (int i = 0; i < n; i += len)
            {
                var w = Complex.One;
                for (int k = 0; k < len / 2; k++)
                {
                    var u = buffer[i + k];
                    var v = buffer[i + k + len / 2] * w;
                    buffer[i + k] = u + v;
                    buffer[i + k + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
    }
}
