using Falcon.Core.Abstractions;

namespace Falcon.Core.Services.Fitting;

internal static class SigmaClip
{
    public static (double[] ClippedData, double LowerBound, double UpperBound) Clip(
        ReadOnlySpan<double> data, double lowSigma = 3.0, double highSigma = 3.0)
    {
        if (data.IsEmpty)
            return (Array.Empty<double>(), double.NegativeInfinity, double.PositiveInfinity);

        var work = data.ToArray();
        double lower, upper;
        while (true)
        {
            double mean = work.Average();
            double std = work.Length > 1
                ? Math.Sqrt(work.Select(v => (v - mean) * (v - mean)).Sum() / (work.Length - 1))
                : 0.0;
            lower = mean - lowSigma * std;
            upper = mean + highSigma * std;
            bool changed = false;
            for (int i = 0; i < work.Length; i++)
            {
                if (work[i] < lower) { work[i] = lower; changed = true; }
                else if (work[i] > upper) { work[i] = upper; changed = true; }
            }
            if (!changed) break;
        }

        var clipped = new List<double>();
        for (int i = 0; i < data.Length; i++)
            if (data[i] > lower && data[i] < upper)
                clipped.Add(data[i]);
        return (clipped.ToArray(), lower, upper);
    }
}

public sealed class CpuChiSquaredEvaluator : IChiSquaredEvaluator
{
    private readonly double[] _chiModels;
    private readonly double[] _chiBuffer;
    private readonly bool[] _sharedClippedMask;
    private readonly double _sigmaClipLow;
    private readonly double _sigmaClipHigh;

    public ReadOnlyMemory<double> ChiModels => _chiModels;
    public int NumModels { get; }
    public int NumWave { get; }
    public bool[] SharedClippedMask => _sharedClippedMask;

    public CpuChiSquaredEvaluator(ReadOnlySpan<double> data, ReadOnlySpan<double> error, ReadOnlySpan<double> models, int numModels, double sigmaClipLow = 3.0, double sigmaClipHigh = 3.0)
    {
        NumModels = numModels;
        NumWave = data.Length;
        _sigmaClipLow = sigmaClipLow;
        _sigmaClipHigh = sigmaClipHigh;
        _chiModels = new double[numModels * NumWave];
        _chiBuffer = new double[NumWave];
        _sharedClippedMask = new bool[NumWave];
        for (int m = 0; m < numModels; m++)
        {
            int row = m * NumWave;
            for (int w = 0; w < NumWave; w++)
                _chiModels[row + w] = (models[row + w] - data[w]) / error[w];
        }
    }

    public void SetSharedClippedMask(ReadOnlySpan<bool> mask)
    {
        mask.CopyTo(_sharedClippedMask);
    }

    public double Evaluate(ReadOnlySpan<double> weights, int branchNumber, Span<bool> clippedMask)
    {
        Array.Clear(_chiBuffer);
        for (int m = 0; m < NumModels; m++)
        {
            if (weights[m] == 0) continue;
            int row = m * NumWave;
            for (int w = 0; w < NumWave; w++)
                _chiBuffer[w] += weights[m] * _chiModels[row + w];
        }

        if (branchNumber == 0)
        {
            var (_, lower, upper) = SigmaClip.Clip(_chiBuffer, _sigmaClipLow, _sigmaClipHigh);
            double sum = 0;
            for (int w = 0; w < NumWave; w++)
            {
                double chi = _chiBuffer[w];
                clippedMask[w] = chi > lower && chi < upper;
                if (clippedMask[w]) sum += chi * chi;
            }
            return sum;
        }

        double chiSq = 0;
        for (int w = 0; w < NumWave; w++)
            if (_sharedClippedMask[w])
            {
                double chi = _chiBuffer[w];
                chiSq += chi * chi;
            }
        return chiSq;
    }
}

public sealed class IterativeSpectralFitter : ISpectralFitter
{
    public FitResult Fit(SpectralFitInput input, FitterOptions? options = null)
    {
        options ??= new FitterOptions();
        input.Validate();
        int numModels = input.NumModels;
        int numWave = input.NumWave;
        double bicTerm = Math.Log(numWave);
        var evaluator = new CpuChiSquaredEvaluator(input.Data.Span, input.Error.Span, input.Models.Span, numModels, options.SigmaClipLow, options.SigmaClipHigh);

        int indexCount = 0;
        var fitList = new List<InternalFit>();
        var zeroWeights = new double[numModels];
        var branch0Mask = new bool[numWave];

        for (int im = 0; im < numModels; im++)
        {
            zeroWeights[im] += 1;
            var weights = (double[])zeroWeights.Clone();
            fitList.Add(CreateFit(evaluator, weights, 0, ref indexCount, branch0Mask, out _));
            zeroWeights[im] -= 1;
        }

        int bestInit = 0;
        double minInitChi = fitList[0].ChiSquared;
        for (int i = 1; i < numModels; i++)
            if (fitList[i].ChiSquared < minInitChi) { minInitChi = fitList[i].ChiSquared; bestInit = i; }

        evaluator.SetSharedClippedMask(fitList[bestInit].ClippedMask ?? branch0Mask);
        int iterateCount = 0;
        bool converged = false;
        fitList = Iterate(fitList, evaluator, numModels, bicTerm, options, ref indexCount, ref iterateCount, out converged);

        fitList.Sort((a, b) => a.ChiSquared.CompareTo(b.ChiSquared));
        double minChi = fitList[0].ChiSquared;
        int mixTop = Math.Min(options.MixTopCount, fitList.Count);
        var extra = Mix(fitList.GetRange(0, mixTop), fitList, evaluator, numModels, options, ref indexCount);
        var total = fitList.Concat(extra).ToList();

        return new FitResult
        {
            Solutions = total.Select(f => new FitSolution
            {
                Index = f.Index,
                BranchNumber = f.BranchNumber,
                Weights = f.Weights,
                ChiSquared = f.ChiSquared,
                ClippedMask = f.ClippedMask
            }).ToList(),
            MinChiSquared = minChi,
            ClippedMask = evaluator.SharedClippedMask,
            IterationCount = iterateCount,
            Converged = converged
        };
    }

    private static List<InternalFit> Iterate(List<InternalFit> fitList, CpuChiSquaredEvaluator evaluator, int numModels, double bicTerm, FitterOptions options, ref int indexCount, ref int iterateCount, out bool converged)
    {
        iterateCount++;
        int countNew = 0;
        int lenList = fitList.Count;
        double previousMinChi = fitList.Min(f => f.ChiSquared);
        for (int f = 0; f < lenList; f++)
        {
            foreach (var child in SpawnChildren(fitList[f], evaluator, numModels, ref indexCount))
            {
                if (child.ChiSquared < previousMinChi - bicTerm)
                {
                    if (options.FitPerIterationCap > 0 && countNew >= options.FitPerIterationCap) break;
                    fitList.Add(child);
                    countNew++;
                }
            }
        }
        if (countNew == 0) { converged = true; return fitList; }
        if (iterateCount >= options.MaxIterations) { converged = false; return fitList; }
        return Iterate(fitList, evaluator, numModels, bicTerm, options, ref indexCount, ref iterateCount, out converged);
    }

    private static List<InternalFit> SpawnChildren(InternalFit parent, CpuChiSquaredEvaluator evaluator, int numModels, ref int indexCount)
    {
        int branchNum = parent.BranchNumber;
        var children = new List<InternalFit>();
        var newWeights = new double[numModels];
        for (int im = 0; im < numModels; im++) newWeights[im] = parent.Weights[im] * branchNum;
        double sumWeights = newWeights.Sum() + 1;
        var scratchMask = new bool[evaluator.NumWave];
        for (int im = 0; im < numModels; im++)
        {
            newWeights[im] += 1;
            var normalized = newWeights.Select(w => w / sumWeights).ToArray();
            children.Add(CreateFit(evaluator, normalized, branchNum, ref indexCount, scratchMask, out _));
            newWeights[im] -= 1;
        }
        return children;
    }

    private static List<InternalFit> Mix(List<InternalFit> topFits, List<InternalFit> fullFitList, CpuChiSquaredEvaluator evaluator, int numModels, FitterOptions options, ref int indexCount)
    {
        var extra = new List<InternalFit>();
        var mixed = new double[numModels];
        var scratchMask = new bool[evaluator.NumWave];
        foreach (var f1 in topFits)
        foreach (var f2 in fullFitList)
        foreach (double q in options.MixQuantities)
        {
            double denom = 1 + q;
            for (int m = 0; m < numModels; m++) mixed[m] = (f1.Weights[m] + q * f2.Weights[m]) / denom;
            extra.Add(CreateFit(evaluator, mixed, f1.BranchNumber + f2.BranchNumber, ref indexCount, scratchMask, out _));
        }
        return extra;
    }

    private static InternalFit CreateFit(CpuChiSquaredEvaluator evaluator, double[] weights, int branchNumber, ref int indexCount, Span<bool> branch0MaskScratch, out bool[]? clippedMask)
    {
        clippedMask = null;
        double chi;
        bool[]? stored = null;
        if (branchNumber == 0)
        {
            chi = evaluator.Evaluate(weights, 0, branch0MaskScratch);
            stored = branch0MaskScratch.ToArray();
            clippedMask = stored;
        }
        else chi = evaluator.Evaluate(weights, branchNumber, Span<bool>.Empty);

        return new InternalFit
        {
            Index = indexCount++,
            BranchNumber = branchNumber,
            Weights = (double[])weights.Clone(),
            ChiSquared = chi,
            ClippedMask = stored
        };
    }

    private sealed class InternalFit
    {
        public required int Index { get; init; }
        public required int BranchNumber { get; init; }
        public required double[] Weights { get; init; }
        public required double ChiSquared { get; init; }
        public bool[]? ClippedMask { get; init; }
    }
}
