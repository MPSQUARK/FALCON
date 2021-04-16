using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MachineLearningSpectralFittingCode
{
    class Spectral_resolution
    {
        // Constructor
        public Spectral_resolution(Accelerator gpu, float[] wave_in, float[] sres_in, bool log10_in = false, string interp_ext_in="extrapolate")
        {
            this.gpu = gpu;

            // wave MUST be 1D and = to sres shape
            if (wave_in.Length != sres_in.Length)
            {
                throw new Exception("wave_in and sres_in must be of equal length");
            }

            this.interpolator = new LinearInterpolation(wave_in, sres_in); // k = 1, and extrapolate 
            this.log10 = log10_in;

            if (log10_in)
            {
                this.dv = this.spectrum_velocity_scale(wave_in);
                return;
            }
            
            // If log10_in is false
            this.dw = wave_in[1] - wave_in[0];

            // min_sig, sig_pd, sig_mask, sig_vo = NONE

        }

        
        // Variable Block
        Accelerator gpu;
        protected bool log10 { get; private set; }
        protected float dv { get; private set; }
        protected float dw { get; private set; }
        /// <summary>
        /// X : wave
        /// Y : sres
        /// </summary>
        public LinearInterpolation interpolator { get; set; }
        protected float min_sig { get; set; }
        public float sig_vo { get; set; }
        public float[] sig_pd { get; set; }
        public bool[] sig_mask { get; set; }


        // Methods
        private float spectrum_velocity_scale(float[] wave)
        {
            return Constants.c_kms * this.spectral_coordinate_step(wave, log:true, _base:MathF.E);
        }

        private float spectral_coordinate_step(float[] wave, bool log=false, float _base=10f)
        {
            if (_base == MathF.E && log)
            {
                return Vector.Diff_Log(gpu, new Vector(wave)).Value.Average();
            }

            if (log)
            {
                return Vector.Diff_LogMult(gpu, new Vector(wave), (1f / MathF.Log(_base))).Value.Average();
            }


            // If log is false
            return Vector.Diff(gpu, new Vector(wave)).Value.Average();
        }

        public void match(Spectral_resolution new_sres, bool no_offset = true, float min_sig_pix = 0f)
        {
            this.GaussianKernelDifference(new_sres,no_offset,min_sig_pix);
        }

        private void GaussianKernelDifference(Spectral_resolution new_sres, bool no_offset = true, float min_sig_pix = 0f)
        {
            this.min_sig = min_sig_pix;

            float[] _wave = this.interpolator.X;
            float[] _sres = this.interpolator.Y;

            float[] interp_sres = new float[_wave.Length];

            // This could be slow??
            for (int i = 0; i < _wave.Length; i++)
            {
                interp_sres[i] = this.interpolator.Interpolate(_wave[i]);
            }

            float[] sig2_vd = DetermineVariance(_wave, interp_sres, _sres);

            // OPTION 1
            if (no_offset)
            {
                this.sig_vo = 0f;

                if (log10)
                {
                    this.finalize_GaussianKernelDifference(Vector.ScalarOperation(gpu, new Vector(sig2_vd), MathF.Pow(this.dv, 2f), "/").Value);
                    return;
                }

                this.finalize_GaussianKernelDifference(Vector.ScalarOperation(gpu, new Vector(sig2_vd), MathF.Pow(this.dw, 2f), "/").Value);
                return;
            }

            // Option 2
            float neg_amin_sig2_vd = -sig2_vd.Min();
            float[] dv = new float[_wave.Length - 1];
            for (int i = 0; i < dv.Length; i++)
            {
                dv[i] = Constants.c_kms * (2 * (_wave[i+1] - _wave[i]) / (_wave[i+1] + _wave[i]) );
            }

            this.sig_vo = neg_amin_sig2_vd - MathF.Pow((this.min_sig * dv.Max()), 2f);

            if (this.sig_vo > 0f)
            {
                sig2_vd = Vector.ScalarOperation(gpu, new Vector(sig2_vd), this.sig_vo, "+").Value;
                this.sig_vo = MathF.Sqrt(this.sig_vo);
            }
            else
            {
                this.sig_vo = 0f;
            }

            float[] sig2_pd = this.convert_vd2pd(sig2_vd, _wave);

            this.finalize_GaussianKernelDifference(sig2_pd);
        }
        private float[] DetermineVariance(float[] wave, float[] interp_sres, float[] sres)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float>(DetermineVarianceKernel);

            MemoryBuffer<float> buffer = gpu.Allocate<float>(wave.Length); // Output
            MemoryBuffer<float> buffer2 = gpu.Allocate<float>(wave.Length); //  Input
            MemoryBuffer<float> buffer3 = gpu.Allocate<float>(interp_sres.Length); //  Input
            MemoryBuffer<float> buffer4 = gpu.Allocate<float>(sres.Length); //  Input

            buffer.MemSetToZero(stream);
            buffer2.MemSetToZero(stream);
            buffer3.MemSetToZero(stream);
            buffer4.MemSetToZero(stream);

            buffer2.CopyFrom(stream, wave, 0, 0, wave.Length);
            buffer3.CopyFrom(stream, interp_sres, 0, 0, interp_sres.Length);
            buffer4.CopyFrom(stream, sres, 0, 0, sres.Length);

            kernelWithStream(stream, buffer.Length, buffer.View, buffer2.View, buffer3.View, buffer4.View, Constants.sig2FWHM, Constants.c_kms);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return Output;
        }
        static void DetermineVarianceKernel(Index1 index, ArrayView<float> Output, ArrayView<float> wave,
            ArrayView<float> interp_sres, ArrayView<float> sres, float fwhm, float c)
        {
            Output[index] = XMath.Pow( (c / wave[index]) ,2f) * XMath.Pow((wave[index] / fwhm), 2f) *
                (1f / XMath.Pow(interp_sres[index], 2f) - 1f / XMath.Pow(sres[index], 2f));
        }

        private float[] convert_vd2pd(float[] sig2_vd, float[] wave)
        {
            if (this.log10)
            {
                return Vector.ScalarOperation(gpu, new Vector(sig2_vd), MathF.Pow(this.dv, -2f), "*").Value;
            }

            float inv_cdw_squ = 1f / MathF.Pow(Constants.c_kms * this.dw,2f);
            float[] Output = new float[sig2_vd.Length];
            for (int i = 0; i < Output.Length; i++)
            {
                Output[i] = sig2_vd[i] * wave[i] * wave[i] * inv_cdw_squ;
            }

            return Output;
        }

        private void finalize_GaussianKernelDifference(float[] sig2_pd)
        {
            int[] indx = UtilityMethods.WhereIsClose(sig2_pd, 0f);
            int[] nindx = UtilityMethods.WhereNot(indx, sig2_pd.Length);
            this.sig_pd = sig2_pd;
            for (int i = 0; i < nindx.Length; i++)
            {
                this.sig_pd[nindx[i]] = sig2_pd[nindx[i]] / MathF.Sqrt(MathF.Abs(sig2_pd[nindx[i]]));
            }
            for (int i = 0; i < indx.Length; i++)
            {
                this.sig_pd[indx[i]] = 0f;
            }

            this.sig_mask = new bool[this.sig_pd.Length];
            for (int i = 0; i < this.sig_mask.Length; i++)
            {
                this.sig_mask[i] = this.sig_pd[i] < -this.min_sig;
            }
            return;
        }

        public float[] adjusted_resolution(int[] indxs )
        {
            float c_by_sig2fwhm = Constants.c_kms / Constants.sig2FWHM;
            float[] output;
            float[] pd2vd;
            if (indxs == null)
            {
                output = new float[this.sig_pd.Length];
                pd2vd = convert_pd2vd(Vector.Power(gpu, new Vector(this.sig_pd), true).Value);

                for (int i = 0; i < output.Length; i++)
                {
                    output[i] = c_by_sig2fwhm / MathF.Sqrt(pd2vd[i])
                         + 1f / MathF.Pow(this.interpolator.Y[i],2f);
                }

                return output;
            }

            output = new float[indxs.Length];
            float[] selected_sig = new float[indxs.Length];
            float[] selected_sres = new float[indxs.Length];
            for (int i = 0; i < output.Length; i++)
            {
                selected_sig[i] = this.sig_pd[indxs[i]] * MathF.Abs(this.sig_pd[indxs[i]]);
                selected_sres[i] = 1f / MathF.Pow(this.interpolator.Y[indxs[i]],2f);
            }

            pd2vd = convert_pd2vd(selected_sig);
            for (int i = 0; i < output.Length; i++)
            {
                output[i] = c_by_sig2fwhm / MathF.Sqrt(pd2vd[i]) + selected_sres[i];
            }

            return output;
        }

        private float[] convert_pd2vd(float[] sig2_pd)
        {
            if (log10)
            {
                return Vector.ScalarOperation(gpu, new Vector(sig2_pd), MathF.Pow(this.dv, 2f), "*").Value;
            }

            float[] sig_sq_cdw = Vector.ScalarOperation(gpu, new Vector(sig2_pd), MathF.Pow(Constants.c_kms * this.dw, 2f), "*").Value;
            float[] output = new float[sig_sq_cdw.Length];
            for (int i = 0; i < output.Length; i++)
            {
                output[i] = sig_sq_cdw[i] / (this.interpolator.X[i] * this.interpolator.X[i]);
            }

            return output;
        }


    }
}
