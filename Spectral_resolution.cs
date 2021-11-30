using BAVCL;
using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FALCON
{
    class Spectral_resolution
    {
        // Constructor
        public Spectral_resolution(GPU gpu, Vector wave_in, Vector sres_in, bool log10_in = false, string interp_ext_in = "extrapolate")
        {
            this.gpu = gpu;

            // wave MUST be 1D and = to sres shape
            if (wave_in.Length != sres_in.Length)
            {
                throw new Exception("wave_in and sres_in must be of equal length");
            }

            this.interpolator = new LinearInterpolation(wave_in.Pull(), sres_in.Pull()); // k = 1, and extrapolate 

            this.log10 = log10_in;

            if (log10_in)
            {
                this.dv = this.spectrum_velocity_scale(wave_in);
                return;
            }

            // If log10_in is false
            this.dw = wave_in.Value[1] - wave_in.Value[0];

            // min_sig, sig_pd, sig_mask, sig_vo = NONE

        }


        // Variable Block
        GPU gpu;
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
        public Vector sig_pd { get; set; }
        public bool[] sig_mask { get; set; }


        // Methods
        private float spectrum_velocity_scale(Vector wave)
        {
            return Constants.c_kms * this.spectral_coordinate_step(wave, log: true, _base: MathF.E);
        }

        private float spectral_coordinate_step(Vector wave, bool log = false, float _base = 10f)
        {
            if (_base == MathF.E && log)
            {
                wave.SyncCPU();
                Vector wave_1 = new Vector(gpu, wave.Value[1..]);
                Vector wave_2 = new Vector(gpu, wave.Value[..^1]);

                return (wave_1.Log_IP(_base) - wave_2.Log_IP(_base)).Mean();
            }

            if (log)
            {
                return Vector.Diff_LogMult(gpu, new Vector(wave), (1f / MathF.Log(_base))).Value.Average();
            }


            // If log is false
            return Vector.Diff(wave).Mean();
        }



        public void match(Spectral_resolution new_sres, bool no_offset = true, float min_sig_pix = 0f)
        {
            this.GaussianKernelDifference(new_sres, no_offset, min_sig_pix);
        }

        private void GaussianKernelDifference(Spectral_resolution new_sres, bool no_offset = true, float min_sig_pix = 0f)
        {
            this.min_sig = min_sig_pix;

            Vector _wave = new Vector(gpu, this.interpolator.X);
            Vector _sres = new Vector(gpu, this.interpolator.Y);

            Vector interp_sres = new Vector(gpu,new float[_wave.Length],1, false);

            // This could be slow??
            for (int i = 0; i < _wave.Length; i++)
            {
                interp_sres.Value[i] = new_sres.interpolator.Interpolate(_wave.Value[i]);
            }
            interp_sres.UpdateCache();

            // Issue?

            //float[] sig2_wd = new float[_wave.Length];
            //for (int i = 0; i < sig2_wd.Length; i++)
            //{
            //    sig2_wd[i] = MathF.Pow(_wave[i] / Constants.sig2FWHM, 2f) * (1.0f / MathF.Pow(interp_sres[i], 2f) - 1.0f / MathF.Pow(_sres[i], 2f));
            //}

            //float[] sig2_vd = new float[sig2_wd.Length];
            //for (int i = 0; i < sig2_vd.Length; i++)
            //{
            //    sig2_vd[i] = MathF.Pow(Constants.c_kms / _wave[i],2f) * sig2_wd[i];
            //}

            Vector sig2_vd = DetermineVariance(_wave, interp_sres, _sres);


            // OPTION 1
            if (no_offset)
            {
                this.sig_vo = 0f;

                if (log10)
                {
                    // ISSUE - dv is Supposed to be 69.02976392336477
                    finalize_GaussianKernelDifference(sig2_vd * XMath.Rcp(this.dv * this.dv));
                    return;
                }

                finalize_GaussianKernelDifference(sig2_vd * XMath.Rcp(this.dw * this.dw));
                return;
            }

            // Option 2

            float neg_amin_sig2_vd = -sig2_vd.Min();
            float[] dv = new float[_wave.Length - 1];
            for (int i = 0; i < dv.Length; i++)
            {
                dv[i] = Constants.c_kms * (2 * (_wave[i + 1] - _wave[i]) / (_wave[i + 1] + _wave[i]));
            }

            this.sig_vo = neg_amin_sig2_vd - MathF.Pow((this.min_sig * dv.Max()), 2f);

            if (this.sig_vo > 0f)
            {
                sig2_vd = sig2_vd + this.sig_vo;
                this.sig_vo = MathF.Sqrt(this.sig_vo);
            }
            else
            {
                this.sig_vo = 0f;
            }


            Vector sig2_pd = this.convert_vd2pd(sig2_vd, _wave);
            finalize_GaussianKernelDifference(sig2_pd);
        }
        private Vector DetermineVariance(Vector wave, Vector interp_sres, Vector sres)
        {
            Vector output = new Vector(gpu, new float[wave.Length]);
            output.IncrementLiveCount();

            wave.IncrementLiveCount();
            interp_sres.IncrementLiveCount();
            sres.IncrementLiveCount();

            var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float>(DetermineVarianceKernel);

            MemoryBuffer<float> outBuffer = output.GetBuffer(); // Output
            MemoryBuffer<float> waveBuffer = wave.GetBuffer(); //  Input
            MemoryBuffer<float> interpsresBuffer = interp_sres.GetBuffer(); //  Input
            MemoryBuffer<float> sresbuffer = sres.GetBuffer(); //  Input

            kernelWithStream(gpu.accelerator.DefaultStream, outBuffer.Length, outBuffer.View, waveBuffer.View, interpsresBuffer.View, sresbuffer.View, Constants.sig2FWHM, Constants.c_kms);

            gpu.accelerator.Synchronize();

            wave.DecrementLiveCount();
            interp_sres.DecrementLiveCount();
            sres.DecrementLiveCount();
            output.DecrementLiveCount();

            return output;
        }
        static void DetermineVarianceKernel(Index1 index, ArrayView<float> Output, ArrayView<float> wave,
            ArrayView<float> interp_sres, ArrayView<float> sres, float fwhm, float c)
        {
            Output[index] = XMath.Pow((c / wave[index]), 2f) * XMath.Pow((wave[index] / fwhm), 2f) *
                (1f / XMath.Pow(interp_sres[index], 2f) - 1f / XMath.Pow(sres[index], 2f));
        }

        private Vector convert_vd2pd(Vector sig2_vd, Vector wave)
        {
            if (this.log10)
            {
                return sig2_vd * MathF.Pow(this.dv, -2f);
            }

            float inv_cdw_squ = XMath.Rcp(Constants.c_kms * Constants.c_kms * this.dw * this.dw);

            return (wave * wave).OP_IP(sig2_vd, Operations.multiply).OP_IP(inv_cdw_squ, Operations.multiply);
        }

        private void finalize_GaussianKernelDifference(Vector sig2_pd)
        {
            int[] indx = UtilityMethods.WhereIsClose(sig2_pd.Pull(), 0f);
            int[] nindx = UtilityMethods.WhereNot(indx, sig2_pd.Length);

            // ISSUE - ALL INF
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

        public float[] adjusted_resolution(int[] indxs)
        {
            
            
            float sig2fwhm_by_c_sq = MathF.Pow(Constants.sig2FWHM / Constants.c_kms, 2f);
            float[] output;
            Vector pd2vd;
            if (indxs.Length == 0)
            {
                output = new float[this.sig_pd.Length];

                pd2vd = convert_pd2vd((+this.sig_pd).OP_IP(this.sig_pd, Operations.multiply));

                for (int i = 0; i < output.Length; i++)
                {
                    output[i] = 1f / MathF.Sqrt(sig2fwhm_by_c_sq * pd2vd.Value[i]
                         + 1f / MathF.Pow(this.interpolator.Y[i], 2f) );
                }

                return output;
            }

            output = new float[indxs.Length];
            Vector selected_sig = new Vector(gpu,new float[indxs.Length],1, false);
            float[] selected_sres = new float[indxs.Length];
            for (int i = 0; i < output.Length; i++)
            {
                selected_sig.Value[i] = this.sig_pd[indxs[i]] * MathF.Abs(this.sig_pd[indxs[i]]);
                selected_sres[i] = 1f / MathF.Pow(this.interpolator.Y[indxs[i]], 2f);
            }
            selected_sig.SyncCPU();

            pd2vd = convert_pd2vd(selected_sig);
            for (int i = 0; i < output.Length; i++)
            {
                output[i] = 1f/ MathF.Sqrt(sig2fwhm_by_c_sq * MathF.Sqrt(pd2vd[i]) + selected_sres[i]);
            }

            return output;
        }

        private Vector convert_pd2vd(Vector sig2_pd)
        {
            if (log10)
            {
                return sig2_pd * (this.dv * this.dv);
            }

            Vector output = sig2_pd * (Constants.c_kms * Constants.c_kms * this.dw * this.dw);
            output.SyncCPU();
            for (int i = 0; i < output.Length; i++)
            {
                output.Value[i] = output.Value[i] * XMath.Rcp(this.interpolator.X[i] * this.interpolator.X[i]);
            }

            return output;
        }


    }

    class VariableGaussianKernel
    {
        public VariableGaussianKernel(GPU gpu, float[] sigma, float minsig=0.01f,int nsig=3, bool integral=false)
        {
            this.gpu = gpu;
            this.n = sigma.Length;
            //this.sigma = (from sig in sigma
            //              select Math.Clamp(sig, minsig, sigma.Max())).ToArray();
            this.sigma = sigma;
            for (int i = 0; i < sigma.Length; i++)
            {
                if (this.sigma[i] < minsig) { this.sigma[i] = minsig; }
                if (this.sigma[i] > sigma.Max()) { this.sigma[i] = sigma.Max(); }
            }

            this.p = (int)MathF.Ceiling(this.sigma.Max() * nsig);
            this.m = 2 * this.p + 1;

            float interval = (2f*MathF.Abs(this.p)) / this.m - 1;
            float[] x2 = (from val in Enumerable.Range(0, this.m)
                                 select MathF.Pow(-this.p + (val * interval), 2f)).ToArray();

            if (!integral)
            {
                this.kernel = kernCalculation(gpu, new Vector(gpu,x2), new Vector(gpu,this.sigma));
            }
            else
            {
                Console.WriteLine("WARNING integral = true in Variable Gaussian KERNEL");
            }

        }

        private int n { get; set; }
        private float[] sigma { get; set; }
        private int p { get; set; }
        private int m { get; set; }
        private Vector kernel { get; set; }
        private GPU gpu;

        public Vector Convolve(GPU gpu, Vector y) //ye=None
        {
            if (y.Length != this.n)
            {
                throw new Exception("Convolution ERROR : y array not the same length as n");
            }

            // assume ye is None
            // if ye = None
            Vector a = this.Create_a(y.Pull());
            Vector result = Vector.MultiplySumAxZero(gpu, a, this.kernel); // PRECISION MAY NEED TO USE DOUBLES

            // else ae = create_a(ye**2) { return sum(a* (this.kernel,axis=0)), sqrt(sum(ae*(this.kernel, axis=0)))}

            return result;
        }

        private Vector Create_a(float[] y)
        {
            List<float> a = new List<float>();

            for (int i = 0; i < m; i++)
            {
                a.AddRange(new float[p]);
                float[] yrange = y[i..(this.n - this.m + i + 1)];
                a.AddRange(yrange);
                a.AddRange(new float[p]);
            }

            return new Vector(gpu, a.ToArray(), this.kernel.Columns);
        }

        private Vector kernCalculation(GPU gpu, Vector x2, Vector sigma)
        {

            x2.IncrementLiveCount();
            sigma.IncrementLiveCount();

            Vector output = new Vector(gpu, new float[x2.Length * sigma.Length], sigma.Length);
            output.IncrementLiveCount();


            var buffer = output.GetBuffer();    // OUTPUT
            var buffer2 = x2.GetBuffer();       // INPUT
            var buffer3 = sigma.GetBuffer();    // INPUT

            var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(kernCalcKERNEL);
            kernelWithStream(gpu.accelerator.DefaultStream, buffer3.Length, buffer.View, buffer2.View, buffer3.View);

            gpu.accelerator.Synchronize();

            x2.DecrementLiveCount();
            sigma.DecrementLiveCount();
            output.DecrementLiveCount();

            return output;
        }

        static void kernCalcKERNEL(Index1 index, ArrayView<float> output, ArrayView<float> x2, ArrayView<float> sig)
        {

            //output[index] = XMath.Exp(x2[XMath.DivRoundDown(index, len)] * 0.5f / XMath.Pow(sig[index % len], 2f));
            //ArrayView<float> partialoutput2 = new ArrayView<float>(x2.Length);
            float sum = 0f;
            int indx = 0;

            for (int i = 0; i < x2.Length; i++)
            {
                indx = index + sig.Length * i;

                output[indx] = XMath.Exp((-x2[i] * 0.5f) / XMath.Pow(sig[index], 2f));
                //partialoutput[i] = XMath.Exp((-x2[i] * 0.5f) / XMath.Pow(sig[index], 2f));
                sum += output[indx];
            }

            sum = 1f / sum; // reciprocal of the sum

            for (int i = 0; i < x2.Length; i++)
            {
                indx = index + sig.Length * i;
                //output[(int)(index + sig.Length * i)] = partialoutput[i] / sum;

                output[indx] *= sum;
            }
            
        }



    }

}
