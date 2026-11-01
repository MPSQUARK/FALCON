using Falcon.Core.Abstractions;
using Falcon.Core.Services.IO;
using ScottPlot;

namespace Falcon.UI.Services;

public static class FalconPlotGenerator
{
    static readonly Color DataBlue = Color.FromHex("#5B9BD5");
    static readonly Color ErrorFill = Color.FromHex("#B8D4EE").WithAlpha(0.35);
    static readonly Color ModelOrange = Color.FromHex("#E0634A");

    const int PlotWidth = 2000;
    const int PlotHeight = 1000;

    public static void SaveSpectrumPlot(
        FalconResultsDocument doc,
        string pngPath,
        int highlightRank = 1,
        PlotRangeMode rangeMode = PlotRangeMode.SnipToModelOverlap)
    {
        string dir = Path.GetDirectoryName(pngPath)!;
        Directory.CreateDirectory(dir);
        string tempPath = Path.Combine(dir, $".{Path.GetFileName(pngPath)}.{Guid.NewGuid():N}.tmp");

        try
        {
            var plt = BuildPlot(doc, highlightRank, rangeMode);
            plt.SavePng(tempPath, PlotWidth, PlotHeight);
            plt.Dispose();

            if (File.Exists(pngPath))
                File.Delete(pngPath);
            File.Move(tempPath, pngPath);
        }
        finally
        {
            if (File.Exists(tempPath))
            {
                try { File.Delete(tempPath); } catch { /* best effort */ }
            }
        }
    }

    static Plot BuildPlot(
        FalconResultsDocument doc,
        int highlightRank,
        PlotRangeMode rangeMode)
    {
        var plt = new Plot();
        plt.FigureBackground.Color = Colors.White;
        plt.DataBackground.Color = Colors.White;

        double[] wave = doc.WavelengthAngstrom;
        double[] data = doc.DataFlux;
        double[] err = doc.ErrorFlux;
        int n = Math.Min(wave.Length, Math.Min(data.Length, err.Length));
        if (n == 0) return plt;

        var selected = doc.TopSolutions.FirstOrDefault(s => s.Rank == highlightRank)
            ?? doc.TopSolutions.FirstOrDefault();

        double[] model = selected is not null && selected.ModelFlux.Length > 0
            ? AlignFlux(wave, selected.ModelFlux)
            : [];

        (int start, int end) = SelectRange(wave, data, model, doc, rangeMode, n);
        wave = wave.AsSpan(start, end - start).ToArray();
        data = data.AsSpan(start, end - start).ToArray();
        err = err.AsSpan(start, end - start).ToArray();
        if (model.Length > 0)
            model = model.AsSpan(start, Math.Min(end - start, model.Length - start)).ToArray();

        var lower = data.Zip(err, (d, e) => d - Math.Abs(e)).ToArray();
        var upper = data.Zip(err, (d, e) => d + Math.Abs(e)).ToArray();
        var band = plt.Add.FillY(wave, lower, upper);
        band.LegendText = "±1σ error";
        band.FillColor = ErrorFill;
        band.LineColor = Colors.Transparent;

        var dataLine = plt.Add.Scatter(wave, data);
        dataLine.LegendText = "Data";
        dataLine.Color = DataBlue;
        dataLine.LineWidth = 0.9f;
        dataLine.MarkerSize = 0;

        if (model.Length > 0)
        {
            var line = plt.Add.Scatter(wave, model);
            line.MarkerSize = 0;
            line.LineWidth = 0.9f;
            line.Color = ModelOrange;
            line.LegendText = $"FALCON (#{selected!.Rank})";
        }

        var (xMin, xMax, yMin, yMax) = ComputeDisplayLimits(wave, data, model);
        plt.Axes.SetLimits(xMin, xMax, yMin, yMax);
        plt.Axes.Margins(0.02, 0.02, 0.02, 0.02);

        double chiDisplay = selected?.ChiSquaredPerDof ?? 0;
        if (chiDisplay <= 0 && doc.DegreesOfFreedom > 0 && doc.MinChiSquared > 0)
            chiDisplay = doc.MinChiSquared / doc.DegreesOfFreedom;

        string baseName = Path.GetFileNameWithoutExtension(doc.SourceFile);
        plt.Title($"χ²/ν: {chiDisplay:F2}  |  {baseName}  |  solution #{selected?.Rank ?? highlightRank}");
        plt.XLabel("Rest-frame wavelength (Å)");
        plt.YLabel("Flux (erg/s/Å/cm²)");
        plt.Grid.MajorLineColor = Color.FromHex("#E8EDF4");
        plt.ShowLegend(Alignment.UpperRight);
        plt.Legend.FontSize = 11;
        plt.Legend.BackgroundColor = Colors.White.WithAlpha(0.92);
        plt.Legend.OutlineColor = Color.FromHex("#D8E0EC");
        return plt;
    }

    static (double XMin, double XMax, double YMin, double YMax) ComputeDisplayLimits(
        double[] wave,
        double[] data,
        double[] model)
    {
        int n = Math.Min(wave.Length, Math.Min(data.Length, model.Length > 0 ? model.Length : data.Length));
        double xMin = wave[0];
        double xMax = wave[Math.Max(0, n - 1)];
        if (n > 1)
        {
            xMin = wave[0];
            xMax = wave[n - 1];
            double xPad = (xMax - xMin) * 0.01;
            xMin -= xPad;
            xMax += xPad;
        }

        if (model.Length == 0)
        {
            var dataOnly = data.Where(double.IsFinite).OrderBy(v => v).ToArray();
            if (dataOnly.Length == 0) return (xMin, xMax, 0, 1);
            double lo = dataOnly[(int)(dataOnly.Length * 0.02)];
            double hi = dataOnly[(int)(dataOnly.Length * 0.98)];
            double pad = Math.Max((hi - lo) * 0.08, 0.1);
            return (xMin, xMax, lo - pad, hi + pad);
        }

        var modelVals = new List<double>(model.Length);
        for (int i = 0; i < model.Length; i++)
            if (double.IsFinite(model[i])) modelVals.Add(model[i]);
        if (modelVals.Count == 0)
            return (xMin, xMax, 0, 1);

        modelVals.Sort();
        double modelLo = modelVals[0];
        double modelHi = modelVals[^1];
        double modelMed = modelVals[modelVals.Count / 2];
        double modelSpan = Math.Max(modelHi - modelLo, Math.Max(modelMed * 0.05, 0.5));
        double band = modelSpan * 3 + Math.Max(modelMed * 0.5, 1);

        var near = new List<double>(n);
        for (int i = 0; i < n; i++)
        {
            if (!double.IsFinite(data[i])) continue;
            double m = i < model.Length && double.IsFinite(model[i]) ? model[i] : modelMed;
            if (Math.Abs(data[i] - m) <= band)
                near.Add(data[i]);
            if (double.IsFinite(m))
                near.Add(m);
        }

        if (near.Count == 0)
            near.AddRange(modelVals);

        near.Sort();
        int count = near.Count;
        double pLo = near[(int)(count * 0.02)];
        double pHi = near[(int)(count * 0.98)];

        pHi = Math.Min(pHi, modelHi + modelSpan * 0.35);
        pLo = Math.Max(pLo, modelLo - modelSpan * 0.15);

        if (pHi <= pLo)
        {
            pLo = modelLo;
            pHi = modelHi;
        }

        double ySpan = pHi - pLo;
        double yPad = ySpan > 0 ? ySpan * 0.06 : Math.Max(Math.Abs(pHi), 1) * 0.06;
        yPad = Math.Max(yPad, modelSpan * 0.08);
        return (xMin, xMax, pLo - yPad, pHi + yPad);
    }

    static (int Start, int End) SelectRange(
        double[] wave,
        double[] data,
        double[] model,
        FalconResultsDocument doc,
        PlotRangeMode rangeMode,
        int n)
    {
        (int dataStart, int dataEnd) = FiniteDataExtent(data, n);

        if (rangeMode == PlotRangeMode.Full)
            return (dataStart, dataEnd);

        int start = dataStart;
        int end = dataEnd;

        if (doc.ModelWaveMin > 0 && doc.ModelWaveMax > doc.ModelWaveMin)
        {
            int modelStart = Array.FindIndex(wave, w => w >= doc.ModelWaveMin);
            int modelEnd = Array.FindLastIndex(wave, w => w <= doc.ModelWaveMax) + 1;
            if (modelStart >= 0 && modelEnd > modelStart)
            {
                start = Math.Max(start, modelStart);
                end = Math.Min(end, modelEnd);
            }
        }
        else if (model.Length > 0)
        {
            double maxAbs = 0;
            for (int i = 0; i < Math.Min(n, model.Length); i++)
                maxAbs = Math.Max(maxAbs, Math.Abs(model[i]));

            double eps = maxAbs * 1e-12;
            int iMin = start;
            int iMax = Math.Min(end, model.Length) - 1;
            while (iMin < iMax && Math.Abs(model[iMin]) <= eps) iMin++;
            while (iMax > iMin && Math.Abs(model[iMax]) <= eps) iMax--;
            start = Math.Max(start, iMin);
            end = Math.Min(end, iMax + 1);
        }

        return start < end ? (start, end) : (dataStart, dataEnd);
    }

    static (int Start, int End) FiniteDataExtent(double[] data, int n)
    {
        int start = 0;
        int end = n;
        for (int i = 0; i < n; i++)
        {
            if (double.IsFinite(data[i]))
            {
                start = i;
                break;
            }
        }

        for (int i = n - 1; i >= start; i--)
        {
            if (double.IsFinite(data[i]))
            {
                end = i + 1;
                break;
            }
        }

        return (start, end);
    }

    static double[] AlignFlux(double[] wave, double[] modelFlux)
    {
        int n = Math.Min(wave.Length, modelFlux.Length);
        var aligned = new double[wave.Length];
        for (int i = 0; i < n; i++)
            aligned[i] = double.IsFinite(modelFlux[i]) ? modelFlux[i] : double.NaN;
        return aligned;
    }
}

/// <summary>Formatted PDF-averaged properties for the results UI.</summary>
public static class ResultsDisplayHelper
{
    public static IReadOnlyList<PdfPropertyGroup> GroupPdfAverages(IReadOnlyDictionary<string, double> av)
    {
        static string Label(string key) => key switch
        {
            "light_age" => "Light-weighted age",
            "light_metal" => "Light-weighted [Z/H]",
            "mass_age" => "Mass-weighted age",
            "mass_metal" => "Mass-weighted [Z/H]",
            "stellar_mass" => "Stellar mass",
            _ => key.Replace('_', ' ')
        };

        static string Units(string key) => key switch
        {
            "light_age" or "mass_age" => "Gyr",
            "light_metal" or "mass_metal" => "dex",
            "stellar_mass" => "G☉",
            _ => ""
        };

        static double DisplayValue(string key, double raw) => key switch
        {
            "light_metal" or "mass_metal" => raw > 0 ? Math.Log10(raw) : double.NaN,
            "stellar_mass" => raw > 0 ? raw / 1e9 : double.NaN,
            _ => raw
        };

        var bases = av.Keys
            .Where(k => !k.Contains("_sig", StringComparison.Ordinal))
            .Where(k => !k.EndsWith("_plus", StringComparison.Ordinal) && !k.EndsWith("_minus", StringComparison.Ordinal))
            .OrderBy(k => k, StringComparer.Ordinal);

        var groups = new List<PdfPropertyGroup>();
        foreach (var key in bases)
        {
            if (!av.TryGetValue(key, out double best)) continue;
            av.TryGetValue($"{key}_1_sig_plus", out double p1);
            av.TryGetValue($"{key}_1_sig_minus", out double m1);
            string? help = key switch
            {
                "stellar_mass" => "10⁹ solar masses (billions of M☉).",
                "light_metal" or "mass_metal" => "log₁₀(metallicity / Z☉); 0 = solar.",
                _ => null
            };
            groups.Add(new PdfPropertyGroup(
                Label(key),
                DisplayValue(key, best),
                Units(key),
                help ?? (p1 > 0 || m1 > 0 ? $"68% interval: +{DisplayValue(key, p1):G4} / −{DisplayValue(key, m1):G4}" : null)));
        }
        return groups;
    }

    public static string FormatSspRatio(IReadOnlyList<FalconSspComponent> components)
    {
        double sum = components.Sum(c => c.LightWeight);
        var parts = components
            .Where(c => c.LightWeight > 1e-6)
            .OrderByDescending(c => c.LightWeight)
            .Take(4)
            .Select(c => $"{c.LightWeight * 100 / (sum > 0 ? sum : 1):F1}% @ {c.AgeGyr:F2} Gyr")
            .ToList();
        return parts.Count == 0 ? "—" : string.Join(" · ", parts);
    }

    public static IReadOnlyList<PdfPropertyGroup> FromSolution(FalconSolutionRecord sol, FalconResultsDocument doc)
    {
        double chiDof = sol.ChiSquaredPerDof;
        if (chiDof <= 0 && doc.DegreesOfFreedom > 0 && doc.MinChiSquared > 0)
            chiDof = doc.MinChiSquared / doc.DegreesOfFreedom;

        double metalDex = sol.LightWeightedMetal > 0 ? Math.Log10(sol.LightWeightedMetal) : double.NaN;
        double stellarMassGmsun = doc.PdfAverages.TryGetValue("stellar_mass", out double mSun) && mSun > 0
            ? mSun / 1e9
            : double.NaN;

        return
        [
            new PdfPropertyGroup("χ² / dof", chiDof, "", "Lower is a better fit to the data."),
            new PdfPropertyGroup("Relative probability", sol.Probability, "", "Higher means this solution is more likely among the explored fits (not an absolute confidence)."),
            new PdfPropertyGroup("Light-weighted age", sol.LightWeightedAgeGyr, "Gyr", null),
            new PdfPropertyGroup("Light-weighted [Z/H]", metalDex, "dex", "log₁₀(metallicity / Z☉); 0 = solar."),
            new PdfPropertyGroup("Stellar mass", stellarMassGmsun, "G☉", "10⁹ solar masses (billions of M☉)."),
        ];
    }

    public static double DisplayMetal(double linearOrLog, bool isLinearLightWeighted = false)
    {
        if (!isLinearLightWeighted) return linearOrLog;
        return linearOrLog > 0 ? Math.Log10(linearOrLog) : double.NaN;
    }
}

public sealed record PdfPropertyGroup(string Label, double Value, string Units, string? HelpText);
