using Falcon.Core.Abstractions;
using Falcon.Core.Numerics;

namespace Falcon.Core.Services.Fitting;

/// <summary>
/// Firefly-equivalent χ²: precomputed float32 <c>chi_models</c> and weighted blend per eval (no per-eval alloc).
/// </summary>
public class FloatChiSquaredEvaluator : IChiSquaredEvaluator
{
    readonly float[] _chiModels;
    readonly float[] _chiBuffer;
    readonly float[] _winsorScratch;
    readonly bool[] _sharedClippedMask;
    readonly double _sigmaClipLow;
    readonly double _sigmaClipHigh;
    readonly double[] _chiModelsDouble;

    public ReadOnlyMemory<double> ChiModels => _chiModelsDouble;
    public int NumModels { get; }
    public int NumWave { get; }
    public bool SupportsBatchEvaluation => true;
    public bool[] SharedClippedMask => _sharedClippedMask;

    public FloatChiSquaredEvaluator(
        ReadOnlySpan<double> data,
        ReadOnlySpan<double> error,
        ReadOnlySpan<double> models,
        int numModels,
        double sigmaClipLow = 3.0,
        double sigmaClipHigh = 3.0)
    {
        NumModels = numModels;
        NumWave = data.Length;
        _sigmaClipLow = sigmaClipLow;
        _sigmaClipHigh = sigmaClipHigh;
        _chiModels = new float[numModels * NumWave];
        _chiModelsDouble = new double[numModels * NumWave];
        _chiBuffer = new float[NumWave];
        _winsorScratch = new float[NumWave];
        _sharedClippedMask = new bool[NumWave];

        var dataF = new float[NumWave];
        var errF = new float[NumWave];
        SimdSpanOps.PrepareDataSpan(data, dataF);
        SimdSpanOps.PrepareErrorSpan(error, errF);
        for (int m = 0; m < numModels; m++)
            SimdSpanOps.BuildChiBankRowMajor(_chiModels, dataF, errF, models, _chiModelsDouble, m, NumWave);
    }

    public FloatChiSquaredEvaluator(
        ReadOnlySpan<float> data,
        ReadOnlySpan<float> error,
        ReadOnlySpan<float> models,
        int numModels,
        double sigmaClipLow = 3.0,
        double sigmaClipHigh = 3.0)
    {
        NumModels = numModels;
        NumWave = data.Length;
        _sigmaClipLow = sigmaClipLow;
        _sigmaClipHigh = sigmaClipHigh;
        _chiModels = new float[numModels * NumWave];
        _chiModelsDouble = new double[numModels * NumWave];
        _chiBuffer = new float[NumWave];
        _winsorScratch = new float[NumWave];
        _sharedClippedMask = new bool[NumWave];

        var errF = new float[NumWave];
        for (int w = 0; w < NumWave; w++)
        {
            float e = error[w];
            errF[w] = !float.IsFinite(e) || e <= 0f ? 1e-30f : e;
        }

        for (int m = 0; m < numModels; m++)
        {
            int row = m * NumWave;
            SimdSpanOps.ComputeChiRow(
                _chiModels.AsSpan(row, NumWave),
                models.Slice(row, NumWave),
                data,
                errF);
            for (int w = 0; w < NumWave; w++)
                _chiModelsDouble[row + w] = _chiModels[row + w];
        }
    }

    public void SetSharedClippedMask(ReadOnlySpan<bool> mask) => mask.CopyTo(_sharedClippedMask);

    public double Evaluate(ReadOnlySpan<double> weights, int branchNumber, Span<bool> clippedMask)
    {
        BlendIntoBuffer(weights, _chiBuffer);
        if (branchNumber == 0)
        {
            if (clippedMask.Length < NumWave)
                throw new ArgumentException("clippedMask must cover all wavelength bins for branch 0.", nameof(clippedMask));
            return ChiSquaredWinsorize.ChiSquaredFromBlended(
                _chiBuffer, _sharedClippedMask, 0, clippedMask, _sigmaClipLow, _sigmaClipHigh);
        }

        double chi = ChiSquaredWinsorize.ChiSquaredFromBlended(
            _chiBuffer, _sharedClippedMask, branchNumber, Span<bool>.Empty, _sigmaClipLow, _sigmaClipHigh);
        if (!double.IsNaN(chi)) return chi;

        var scratch = new bool[NumWave];
        return Evaluate(weights, 0, scratch);
    }

    public int EvaluateBatch(
        ReadOnlySpan<double> flatWeights,
        int batchSize,
        ReadOnlySpan<int> branchNumbers,
        Span<double> chiSquaredOut)
    {
        var scratchMask = new bool[NumWave];
        for (int b = 0; b < batchSize; b++)
        {
            ReadOnlySpan<double> weights = flatWeights.Slice(b * NumModels, NumModels);
            chiSquaredOut[b] = Evaluate(weights, branchNumbers[b], branchNumbers[b] == 0 ? scratchMask : Span<bool>.Empty);
        }
        return batchSize;
    }

    void BlendIntoBuffer(ReadOnlySpan<double> weights, Span<float> dest)
    {
        dest.Clear();
        for (int m = 0; m < NumModels; m++)
        {
            if (weights[m] == 0) continue;
            SimdSpanOps.AddWeighted(dest, _chiModels.AsSpan(m * NumWave, NumWave), (float)weights[m]);
        }
    }

    public void Dispose() { }
}

/// <summary>Backward-compatible alias for tests and callers.</summary>
public sealed class CpuChiSquaredEvaluator : FloatChiSquaredEvaluator
{
    public CpuChiSquaredEvaluator(
        ReadOnlySpan<double> data,
        ReadOnlySpan<double> error,
        ReadOnlySpan<double> models,
        int numModels,
        double sigmaClipLow = 3.0,
        double sigmaClipHigh = 3.0)
        : base(data, error, models, numModels, sigmaClipLow, sigmaClipHigh)
    {
    }
}
