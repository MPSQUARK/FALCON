using BAVCL;
using BAVCL.Core;
using Falcon.Core.Abstractions;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Numerics;
using ILGPU;
using ILGPU.Runtime;

namespace Falcon.Core.Services.Fitting;

/// <summary>
/// GPU χ² blend: persistent 2D chi bank (wave × model), reused weight/batch buffers, batched kernel for Mix.
/// Sigma-clipping runs on CPU after blend (firefly <c>sigmaclip</c> parity).
/// </summary>
public sealed class GpuChiSquaredEvaluator : IChiSquaredEvaluator
{
    const int DefaultBatchCapacity = 128;

    readonly GPU _gpu;
    readonly Vector _chiBank;
    readonly Vector _weightVec;
    readonly Vector _weightBatchVec;
    readonly Vector _blendedBatch;
    readonly float[] _weightScratch;
    readonly float[] _weightBatchScratch;
    readonly float[] _chiBuffer;
    readonly float[] _batchBlendedCpu;
    readonly bool[] _sharedClippedMask;
    readonly double[] _chiModelsDouble;
    readonly double _sigmaClipLow;
    readonly double _sigmaClipHigh;
    readonly int _batchCapacity;

    public ReadOnlyMemory<double> ChiModels => _chiModelsDouble;
    public int NumModels { get; }
    public int NumWave { get; }
    public bool SupportsBatchEvaluation => true;
    public bool[] SharedClippedMask => _sharedClippedMask;

    public GpuChiSquaredEvaluator(
        GPU gpu,
        ReadOnlySpan<double> data,
        ReadOnlySpan<double> error,
        ReadOnlySpan<double> models,
        int numModels,
        double sigmaClipLow = 3.0,
        double sigmaClipHigh = 3.0,
        int batchCapacity = DefaultBatchCapacity)
    {
        FalconGpuHost.EnsureReady(gpu);
        _gpu = gpu;
        NumModels = numModels;
        NumWave = data.Length;
        _sigmaClipLow = sigmaClipLow;
        _sigmaClipHigh = sigmaClipHigh;
        _batchCapacity = Math.Max(1, batchCapacity);
        _chiModelsDouble = new double[numModels * NumWave];
        _chiBuffer = new float[NumWave];
        _weightScratch = new float[numModels];
        _weightBatchScratch = new float[_batchCapacity * numModels];
        _batchBlendedCpu = new float[_batchCapacity * NumWave];
        _sharedClippedMask = new bool[NumWave];

        var bank = new float[NumWave * numModels];
        var dataF = new float[NumWave];
        var errF = new float[NumWave];
        var modelF = new float[NumWave];
        var chiRow = new float[NumWave];
        SimdSpanOps.PrepareDataSpan(data, dataF);
        SimdSpanOps.PrepareErrorSpan(error, errF);
        for (int m = 0; m < numModels; m++)
        {
            int row = m * NumWave;
            ReadOnlySpan<double> modelRow = models.Slice(row, NumWave);
            for (int w = 0; w < NumWave; w++)
            {
                double mod = modelRow[w];
                modelF[w] = double.IsFinite(mod) ? (float)mod : 0f;
            }

            SimdSpanOps.ComputeChiRow(chiRow, modelF, dataF, errF);
            SimdSpanOps.ScatterChiRowToGpuLayout(chiRow, bank, m, numModels, NumWave);
            for (int w = 0; w < NumWave; w++)
                _chiModelsDouble[row + w] = chiRow[w];
        }

        _chiBank = new Vector(gpu, bank, numModels);
        _weightVec = new Vector(gpu, _weightScratch);
        _weightBatchVec = new Vector(gpu, _weightBatchScratch);
        _blendedBatch = new Vector(gpu, new float[_batchCapacity * NumWave]);
    }

    public void SetSharedClippedMask(ReadOnlySpan<bool> mask) => mask.CopyTo(_sharedClippedMask);

    public double Evaluate(ReadOnlySpan<double> weights, int branchNumber, Span<bool> clippedMask)
    {
        for (int m = 0; m < NumModels; m++)
            _weightScratch[m] = (float)weights[m];
        _weightVec.UpdateCache(_weightScratch);
        LaunchBlend(_weightVec.GetBuffer().View, 1);
        PullBlendedRow(0, _chiBuffer);
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
        if (batchSize <= 0) return 0;

        FalconGpuHost.EnsureReady(_gpu);
        var scratchMask = new bool[NumWave];

        for (int offset = 0; offset < batchSize; offset += _batchCapacity)
        {
            int chunk = Math.Min(_batchCapacity, batchSize - offset);
            for (int b = 0; b < chunk; b++)
            {
                int baseW = b * NumModels;
                int src = (offset + b) * NumModels;
                for (int m = 0; m < NumModels; m++)
                    _weightBatchScratch[baseW + m] = (float)flatWeights[src + m];
            }

            _weightBatchVec.UpdateCache(_weightBatchScratch);
            LaunchBlend(_weightBatchVec.GetBuffer().View, chunk);
            PullBlendedRows(chunk);

            for (int b = 0; b < chunk; b++)
            {
                int global = offset + b;
                ReadOnlySpan<float> row = _batchBlendedCpu.AsSpan(b * NumWave, NumWave);
                ReadOnlySpan<double> weights = flatWeights.Slice(global * NumModels, NumModels);
                int branch = branchNumbers[global];
                if (branch == 0)
                    chiSquaredOut[global] = Evaluate(weights, 0, scratchMask);
                else
                    chiSquaredOut[global] = Evaluate(weights, branch, Span<bool>.Empty);
            }
        }

        return batchSize;
    }

    void LaunchBlend(ArrayView<float> weightRows, int batchSize)
    {
        using (GpuScope.Begin(_blendedBatch, _chiBank))
        {
            ChiBlendGpu.BlendBatch(
                _gpu,
                _chiBank.GetBuffer().View,
                weightRows,
                _blendedBatch.GetBuffer().View,
                batchSize,
                NumModels,
                NumWave);
        }
    }

    void PullBlendedRow(int batchIndex, Span<float> dest)
    {
        using (GpuScope.BeginReadOnly(_blendedBatch))
            _blendedBatch.RetrieveReadOnlySpan().Slice(batchIndex * NumWave, NumWave).CopyTo(dest);
    }

    void PullBlendedRows(int chunk)
    {
        using (GpuScope.BeginReadOnly(_blendedBatch))
            _blendedBatch.RetrieveReadOnlySpan().Slice(0, chunk * NumWave).CopyTo(_batchBlendedCpu.AsSpan(0, chunk * NumWave));
    }

    public void Dispose() { }
}
