# BAVCL Feature Requests (from FALCON)

Track APIs that would let FALCON delegate more work to BAVCL.

| Priority | Feature | FALCON usage today |
|----------|---------|-------------------|
| High | GPU FFT / IFFT | HPF dust determination (`DustService` uses CPU FFT) |
| High | GPU convolution (variable Gaussian kernel) | Spectral resolution downgrade |
| Medium | Sigma-clip reduction on GPU | Iterative fitter residual clipping |
| Medium | GPU median / percentile | `NormaliseSpec` hot path |
| Medium | GPU optimisation (PSO / grid search) | E(B-V) grid search |
| Low | FITS I/O | CSharpFITS + Python scripts |
| Low | HDF5 I/O | HDF5-CSharp + Python scripts |

## Notes

- FALCON loads `KernelWorkloads.Statistics` explicitly at startup.
- `InternalsVisibleTo Falcon.Core` added to BAVCL for transitional `Vector.Value` access during migration.
