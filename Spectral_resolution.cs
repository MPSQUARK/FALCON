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

            // call this.interpolator = InterpolatedUnivariateSpline(wave,sres,k=1,ext=interp_ext)
            this.log10 = log10_in;
            // cnst = Constants()
            // c = Constants.c_kms


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

    }
}
