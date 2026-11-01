using Falcon.Core.Abstractions;
using Falcon.Core.Diagnostics;
using BAVCL;

namespace Falcon.Core.Services.Fitting;

public sealed class IterativeSpectralFitter(GPU? gpu = null) : ISpectralFitter
{
    static readonly double[] DefaultMixQuantities = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0];

    public FitResult Fit(SpectralFitInput input, FitterOptions? options = null)
    {
        options ??= new FitterOptions();
        input.Validate();
        int numModels = input.NumModels;
        int numWave = input.NumWave;
        double bicTerm = Math.Log(numWave);

        using IChiSquaredEvaluator evaluator = CreateEvaluator(input, numModels, options);

        return FitCore(evaluator, numModels, numWave, bicTerm, options);
    }

    IChiSquaredEvaluator CreateEvaluator(SpectralFitInput input, int numModels, FitterOptions options)
    {
        if (options.UseGpuEvaluator && gpu?.accelerator is not null)
        {
            return new GpuChiSquaredEvaluator(
                gpu,
                input.Data.Span,
                input.Error.Span,
                input.Models.Span,
                numModels,
                options.SigmaClipLow,
                options.SigmaClipHigh);
        }

        return new FloatChiSquaredEvaluator(
            input.Data.Span,
            input.Error.Span,
            input.Models.Span,
            numModels,
            options.SigmaClipLow,
            options.SigmaClipHigh);
    }

    static FitResult FitCore(
        IChiSquaredEvaluator evaluator,
        int numModels,
        int numWave,
        double bicTerm,
        FitterOptions options)
    {
        int indexCount = 0;
        int mixEvalCount = 0;
        var fitList = new List<InternalFit>();
        var seenRatios = new HashSet<string>();
        var zeroWeights = new double[numModels];
        var branch0Mask = new bool[numWave];

        for (int im = 0; im < numModels; im++)
        {
            zeroWeights[im] += 1;
            var weights = (double[])zeroWeights.Clone();
            var fit = TryCreateFit(evaluator, weights, 0, ref indexCount, branch0Mask, seenRatios);
            if (fit is not null) fitList.Add(fit);
            zeroWeights[im] -= 1;
        }

        int bestInit = 0;
        double minInitChi = fitList[0].ChiSquared;
        for (int i = 1; i < fitList.Count; i++)
            if (fitList[i].ChiSquared < minInitChi) { minInitChi = fitList[i].ChiSquared; bestInit = i; }

        evaluator.SetSharedClippedMask(fitList[bestInit].ClippedMask ?? branch0Mask);
        int iterateCount = 0;
        bool converged = false;
        fitList = Iterate(fitList, evaluator, numModels, bicTerm, options, ref indexCount, ref iterateCount, seenRatios, out converged);

        fitList.Sort((a, b) => a.ChiSquared.CompareTo(b.ChiSquared));
        double minChi = fitList[0].ChiSquared;
        int mixTop = Math.Min(options.MixTopCount, fitList.Count);
        var extra = Mix(fitList.GetRange(0, mixTop), fitList, evaluator, numModels, options, ref indexCount, ref mixEvalCount, seenRatios);
        var total = fitList.Concat(extra).OrderBy(f => f.ChiSquared).ToList();

        EmitFitSummary(indexCount, iterateCount, mixEvalCount, fitList.Count, converged);

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

    static void EmitFitSummary(int indexCount, int iterateCount, int mixEvalCount, int finalListSize, bool converged)
    {
        if (!FitTrace.Enabled) return;
        FitTrace.Emit("fitter_summary", new Dictionary<string, object>
        {
            ["index_count"] = indexCount,
            ["iterate_count"] = iterateCount,
            ["mix_eval_count"] = mixEvalCount,
            ["final_fit_list_size"] = finalListSize,
            ["converged"] = converged,
        });
    }

    private static List<InternalFit> Iterate(
        List<InternalFit> fitList,
        IChiSquaredEvaluator evaluator,
        int numModels,
        double bicTerm,
        FitterOptions options,
        ref int indexCount,
        ref int iterateCount,
        HashSet<string> seenRatios,
        out bool converged)
    {
        iterateCount++;
        int countNew = 0;
        int lenList = fitList.Count;
        double previousMinChi = fitList.Min(f => f.ChiSquared);
        for (int f = 0; f < lenList; f++)
        {
            foreach (var child in SpawnChildren(fitList[f], evaluator, numModels, iterateCount, ref indexCount, seenRatios, options.MaxSspComponents))
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
        return Iterate(fitList, evaluator, numModels, bicTerm, options, ref indexCount, ref iterateCount, seenRatios, out converged);
    }

    private static List<InternalFit> SpawnChildren(
        InternalFit parent,
        IChiSquaredEvaluator evaluator,
        int numModels,
        int iterationBranch,
        ref int indexCount,
        HashSet<string> seenRatios,
        int maxSspComponents)
    {
        var children = new List<InternalFit>();
        var newWeights = new double[numModels];
        for (int im = 0; im < numModels; im++) newWeights[im] = parent.Weights[im] * iterationBranch;
        double sumWeights = newWeights.Sum() + 1;
        var scratchMask = new bool[evaluator.NumWave];
        int activeParent = FitSolutionFingerprint.CountActiveComponents(parent.Weights);
        for (int im = 0; im < numModels; im++)
        {
            if (parent.Weights[im] <= 1e-12 && activeParent >= maxSspComponents)
                continue;

            newWeights[im] += 1;
            var normalized = newWeights.Select(w => w / sumWeights).ToArray();
            var child = TryCreateFit(evaluator, normalized, iterationBranch, ref indexCount, scratchMask, seenRatios);
            if (child is not null) children.Add(child);
            newWeights[im] -= 1;
        }
        return children;
    }

    private static List<InternalFit> Mix(
        List<InternalFit> topFits,
        List<InternalFit> fullFitList,
        IChiSquaredEvaluator evaluator,
        int numModels,
        FitterOptions options,
        ref int indexCount,
        ref int mixEvalCount,
        HashSet<string> seenRatios)
    {
        if (options.ExperimentalMixSearch)
            return MixGoldenSection(topFits, fullFitList, evaluator, numModels, options, ref indexCount, ref mixEvalCount, seenRatios);

        return MixFireflyExact(topFits, fullFitList, evaluator, numModels, options, ref indexCount, ref mixEvalCount, seenRatios);
    }

    private static List<InternalFit> MixFireflyExact(
        List<InternalFit> topFits,
        List<InternalFit> fullFitList,
        IChiSquaredEvaluator evaluator,
        int numModels,
        FitterOptions options,
        ref int indexCount,
        ref int mixEvalCount,
        HashSet<string> seenRatios)
    {
        var candidates = new List<(double[] Weights, int Branch)>();
        var mixed = new double[numModels];
        foreach (var f1 in topFits)
        foreach (var f2 in fullFitList)
        foreach (double q in options.MixQuantities)
        {
            double denom = 1 + q;
            for (int m = 0; m < numModels; m++) mixed[m] = (f1.Weights[m] + q * f2.Weights[m]) / denom;
            if (FitSolutionFingerprint.CountActiveComponents(mixed) > options.MaxSspComponents)
                continue;
            string ratioKey = FitSolutionFingerprint.FromWeights(mixed);
            if (string.IsNullOrEmpty(ratioKey) || !seenRatios.Add(ratioKey))
                continue;
            candidates.Add(((double[])mixed.Clone(), f1.BranchNumber + f2.BranchNumber));
        }

        return ScoreMixCandidates(evaluator, candidates, ref indexCount, ref mixEvalCount);
    }

    private static List<InternalFit> MixGoldenSection(
        List<InternalFit> topFits,
        List<InternalFit> fullFitList,
        IChiSquaredEvaluator evaluator,
        int numModels,
        FitterOptions options,
        ref int indexCount,
        ref int mixEvalCount,
        HashSet<string> seenRatios)
    {
        var extra = new List<InternalFit>();
        var mixed = new double[numModels];
        const double qMin = 1e-7;
        const double qMax = 1.0;
        const double phi = 0.6180339887498949;
        int maxEvals = Math.Max(8, options.ExperimentalMixSearchMaxEvals);

        foreach (var f1 in topFits)
        foreach (var f2 in fullFitList)
        {
            if (!IsTwoComponentPair(f1.Weights, f2.Weights, numModels))
            {
                foreach (double q in DefaultMixQuantities)
                {
                    BlendWeights(f1.Weights, f2.Weights, q, mixed);
                    if (FitSolutionFingerprint.CountActiveComponents(mixed) > options.MaxSspComponents) continue;
                    var fit = TryCreateFit(evaluator, mixed, f1.BranchNumber + f2.BranchNumber, ref indexCount, new bool[evaluator.NumWave], seenRatios);
                    if (fit is not null) { extra.Add(fit); mixEvalCount++; }
                }
                continue;
            }

            double lo = qMin;
            double hi = qMax;
            int evals = 0;
            while (evals < maxEvals && hi - lo > 1e-6)
            {
                double q1 = hi - phi * (hi - lo);
                double q2 = lo + phi * (hi - lo);
                BlendWeights(f1.Weights, f2.Weights, q1, mixed);
                var fit1 = TryCreateFit(evaluator, mixed, f1.BranchNumber + f2.BranchNumber, ref indexCount, new bool[evaluator.NumWave], seenRatios);
                mixEvalCount++;
                evals++;
                BlendWeights(f1.Weights, f2.Weights, q2, mixed);
                var fit2 = TryCreateFit(evaluator, mixed, f1.BranchNumber + f2.BranchNumber, ref indexCount, new bool[evaluator.NumWave], seenRatios);
                mixEvalCount++;
                evals++;
                if (fit1 is null && fit2 is null) break;
                double chi1 = fit1?.ChiSquared ?? double.MaxValue;
                double chi2 = fit2?.ChiSquared ?? double.MaxValue;
                if (chi1 < chi2)
                {
                    hi = q2;
                    if (fit1 is not null) extra.Add(fit1);
                }
                else
                {
                    lo = q1;
                    if (fit2 is not null) extra.Add(fit2);
                }
            }
        }

        return extra;
    }

    static bool IsTwoComponentPair(double[] w1, double[] w2, int numModels)
    {
        int active = 0;
        for (int m = 0; m < numModels; m++)
            if (w1[m] > 1e-12 || w2[m] > 1e-12) active++;
        return active == 2;
    }

    static void BlendWeights(double[] w1, double[] w2, double q, Span<double> dest)
    {
        double denom = 1 + q;
        for (int m = 0; m < dest.Length; m++)
            dest[m] = (w1[m] + q * w2[m]) / denom;
    }

    private static List<InternalFit> ScoreMixCandidates(
        IChiSquaredEvaluator evaluator,
        List<(double[] Weights, int Branch)> candidates,
        ref int indexCount,
        ref int mixEvalCount)
    {
        var extra = new List<InternalFit>();
        if (candidates.Count == 0) return extra;

        if (!evaluator.SupportsBatchEvaluation || candidates.Count < 4)
        {
            var singleMask = new bool[evaluator.NumWave];
            foreach (var (weights, branch) in candidates)
            {
                var fit = TryCreateFitAlreadyDeduped(evaluator, weights, branch, ref indexCount, singleMask);
                mixEvalCount++;
                if (fit is not null) extra.Add(fit);
            }
            return extra;
        }

        int n = candidates.Count;
        var flat = new double[n * evaluator.NumModels];
        var branches = new int[n];
        var chi = new double[n];
        for (int i = 0; i < n; i++)
        {
            candidates[i].Weights.CopyTo(flat.AsSpan(i * evaluator.NumModels, evaluator.NumModels));
            branches[i] = candidates[i].Branch;
        }

        evaluator.EvaluateBatch(flat, n, branches, chi);
        mixEvalCount += n;

        var scratchMask = new bool[evaluator.NumWave];
        for (int i = 0; i < n; i++)
        {
            bool[]? stored = null;
            if (candidates[i].Branch == 0)
            {
                evaluator.Evaluate(candidates[i].Weights, 0, scratchMask);
                stored = scratchMask.ToArray();
            }

            extra.Add(new InternalFit
            {
                Index = indexCount++,
                BranchNumber = candidates[i].Branch,
                Weights = candidates[i].Weights,
                ChiSquared = chi[i],
                ClippedMask = stored
            });
        }

        return extra;
    }

    private static InternalFit? TryCreateFit(
        IChiSquaredEvaluator evaluator,
        double[] weights,
        int branchNumber,
        ref int indexCount,
        Span<bool> branch0MaskScratch,
        HashSet<string> seenRatios)
    {
        string ratioKey = FitSolutionFingerprint.FromWeights(weights);
        if (string.IsNullOrEmpty(ratioKey) || !seenRatios.Add(ratioKey))
            return null;
        return TryCreateFitAlreadyDeduped(evaluator, weights, branchNumber, ref indexCount, branch0MaskScratch);
    }

    private static InternalFit? TryCreateFitAlreadyDeduped(
        IChiSquaredEvaluator evaluator,
        double[] weights,
        int branchNumber,
        ref int indexCount,
        Span<bool> branch0MaskScratch)
    {
        bool[]? clippedMask = null;
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
