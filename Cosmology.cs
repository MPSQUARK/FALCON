using ILGPU;
using ILGPU.Runtime;
using ILGPU.Algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    class Cosmology
    {

        float
            Oc0,
            Ob0,
            Om0,
            Ode0,
            H0,
            n,
            sigma8,
            tau,
            z_reion,
            t0,
            Tcmb0,
            Neff;
        bool flat = true;
        float[] m_nu;
        double
            H0_s,
            critical_density0, 
            Ogamma0;


        float
            Odm0,
            Onu0,
            Ok0,
            Tnu0,
            h,
            hubble_distance,
            hubble_time,
            neff_per_nu;
        int
            nneutrinos,
            nmasslessnu,
            nmassivenu;
        bool massivenu;
        float[]
            massivenu_mass,
            nu_y;

        private Func<float, inv_efunc_scalar_args_struct, float> inv_efunc_scalar;
        private inv_efunc_scalar_args_struct inv_efunc_scalar_args;

        public Cosmology()
        {
            this.Oc0 = 0.2589f;
            this.Ob0 = 0.04860f;
            this.Om0 = 0.3075f;
            this.H0 = 67.74f;
            this.n = 0.9667f;
            this.sigma8 = 0.8159f;
            this.tau = 0.066f;
            this.z_reion = 8.8f;
            this.t0 = 13.799f;
            this.Tcmb0 = 2.7255f;
            this.Neff = 3.046f;
            this.flat = true;
            this.m_nu = new float[] { 0f, 0f, 0.06f };

        }

        public void Initialise() 
        {
            this.Odm0 = this.Om0 - this.Ob0;
            this.h = this.H0 * 0.01f;
            this.hubble_distance = (Constants.c * 0.001f) / this.H0;
            this.H0_s = this.H0 * Constants.H0units_to_invs;
            this.hubble_time = (float)(Constants.sec_to_Gyr / this.H0_s);
            this.critical_density0 = Constants.critdens_const * Math.Pow(this.H0_s, 2);
            this.nneutrinos = (int)MathF.Floor(this.Neff);

            this.massivenu = false;
            if (this.nneutrinos > 0 && this.Tcmb0 > 0) 
            {
                this.neff_per_nu = this.Neff / this.nneutrinos;

                if (this.m_nu.Max() == 0)
                {
                    this.nmasslessnu = this.nneutrinos;
                    this.nmassivenu = 0;
                }
                else
                {
                    this.massivenu = true;
                    if (this.m_nu.Count() != this.nneutrinos)
                    {
                        throw new Exception("Class Cosmology : func Initialise - Unexpected number of neutrino masses");
                    }

                    this.massivenu_mass = (from nu in this.m_nu
                                           where nu != 0
                                           select nu).ToArray();
                    this.nmassivenu = this.massivenu_mass.Count();
                    this.nmasslessnu = this.nneutrinos - this.nmassivenu;
                }
            }

            if (this.Tcmb0 > 0)
            {
                this.Ogamma0 = (Constants.a_B_c2 * Math.Pow(this.Tcmb0, 4d) / this.critical_density0);
                this.Tnu0 = 0.7137658555036082f * this.Tcmb0;

                if (this.massivenu)
                {
                    int length = this.massivenu_mass.Count();
                    this.nu_y = new float[length];

                    for (int i = 0; i < length; i++)
                    {
                        this.nu_y[i] = (float)(this.massivenu_mass[i] / (Constants.kB_evK * this.Tnu0));
                    }


                    this.Onu0 = (float)(this.Ogamma0 * nu_relative_density(0f));

                }
                else
                {
                    this.Onu0 = (float)(0.22710731766f * this.Neff * this.Ogamma0);
                }
            }
            else
            {
                throw new Exception("Class Cosmology : func Initialise - Unexpected value of Tcmb0, should be tcmb0 > 0 Kelvin");
            }

            this.Ode0 = (float)(1f - (this.Om0 + this.Ogamma0 + this.Onu0));
            this.Ok0 = 0; // Assuming no curvature

            if (!this.massivenu)
            {
                this.inv_efunc_scalar = flcdm_inv_efunc_nomnu;
                this.inv_efunc_scalar_args = new inv_efunc_scalar_args_struct
                {
                    Om0 = this.Om0,
                    Ode0 = this.Ode0,
                    Or0 = (float)(this.Ogamma0 + this.Onu0)
                };
            }
            else
            {
                this.inv_efunc_scalar = flcdm_inv_efunc;
                this.inv_efunc_scalar_args = new inv_efunc_scalar_args_struct
                {
                    Om0 = this.Om0,
                    Ode0 = this.Ode0,
                    Ogamma0 = (float)this.Ogamma0,
                    neff_per_nu = this.neff_per_nu,
                    nmasslessnu = this.nmasslessnu,
                    nu_y = this.nu_y
                };
            }

            // Further optimisations




        }


        private struct inv_efunc_scalar_args_struct
        {
            public float Om0;
            public float Ode0;
            public float Or0;

            public float Ogamma0;
            public float neff_per_nu;
            public float nmasslessnu;
            public float[] nu_y;
        }


        private float nu_relative_density(float z)
        {
            // See Komatsu et al. 2011, eq 26 and the surrounding discussion
            float prefac = 0.22710731766f;  // 7/8 (4/11)^4/3 -- see any cosmo book

            // Assume that z is always a scalar
            if (!this.massivenu)
            {
                return prefac * this.Neff;
            }

            // These are purely fitting constants -- see the Komatsu paper
            float p = 1.83f;
            float invp = 0.54644808743f;  // 1.0 / p
            float k = 0.3173f;

            float[] rel_mass_per = new float[this.nu_y.Count()];
            for (int i = 0; i < this.nu_y.Count(); i++)
            {
                rel_mass_per[i] = MathF.Pow(1f + MathF.Pow(k * (this.nu_y[i] / (1f + z)), p), invp);
            }
            float rel_mass = rel_mass_per.Sum() + this.nmasslessnu;

            return prefac * this.neff_per_nu * rel_mass;
        }

        private float flcdm_inv_efunc_nomnu(float z, inv_efunc_scalar_args_struct args)
        {
            float opz = 1 + z;
            return MathF.Pow( MathF.Pow(opz,3f) * (opz * args.Or0 + args.Om0) + args.Ode0, -0.5f);
        }

        private float flcdm_inv_efunc(float z, inv_efunc_scalar_args_struct args)
        {
            float opz = 1f + z;
            float Or0 = args.Ogamma0 * (1f + nufunc(opz, args.neff_per_nu, args.nmasslessnu, args.nu_y));
            return MathF.Pow( MathF.Pow(opz,3f) * (opz * Or0 + args.Om0) + args.Ode0, -0.5f);
        }

        private float nufunc(float opz, float neff_per_nu, float nmasslessnu, float[] nu_y)
        {
            int N = nu_y.Count();
            float k = 0.3173f / opz;
            float rel_mass_sum = nmasslessnu;

            if (N == 1)
            {
                return 0.22710731766f * neff_per_nu * (rel_mass_sum + MathF.Pow(1f + MathF.Pow(k * nu_y[0], 1.83f), 0.54644808743f));
            }
            else
            {
                for (int i = 0; i < N; i++)
                {
                    rel_mass_sum += MathF.Pow(1f + MathF.Pow(k * nu_y[i], 1.83f), 0.54644808743f);
                }
                return 0.22710731766f * neff_per_nu * rel_mass_sum;
            }

        }
 
        public float luminosity_distance(AcceleratorId acceleratorId, float redshift)
        {
            return (1f + redshift) * comoving_transverse_distance(acceleratorId, redshift);
        }

        public float comoving_transverse_distance(AcceleratorId acceleratorId, float redshift)
        {
            float dc = integral_comoving_distance(acceleratorId, redshift);
            return dc;
        }

        public float integral_comoving_distance(AcceleratorId acceleratorId, float redshift)
        {
            return this.hubble_distance * GPU_Integration(acceleratorId, redshift, 1e-8f, this.inv_efunc_scalar_args); //Integrate(inv_efunc_scalar, redshift, 1e-8f, inv_efunc_scalar_args);
        }


        private float Integrate(Func<float, inv_efunc_scalar_args_struct, float> func, float z, float da, inv_efunc_scalar_args_struct args )
        {
            int itter = (int)(z / da);
            float[] vals = new float[itter];

            for (int i = 0; i < itter; i++)
            {
                vals[i] = func(i*da, args) * da;
            }

            return vals.Sum();
        }


        // Integrates the cosmology func to get luminosity distance
        private float GPU_Integration(AcceleratorId acceleratorId, float z, float dz, inv_efunc_scalar_args_struct args)
        {
            int length = (int)(z / dz);

            using var context = new Context();
            context.EnableAlgorithms();

            using var accelerator = Accelerator.Create(context, acceleratorId);

            

            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView< float >, float, float, float, float, float, float, float>(GPU_IntegrationKernal);

            #region


            var buffer = accelerator.Allocate<float>(length);
            buffer.MemSetToZero();

            kernel(buffer.Length, buffer.View, dz, args.Ogamma0, args.Om0, args.Ode0, args.neff_per_nu, args.nmasslessnu, args.nu_y[0]);

            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();

            return Output.Sum();

            #endregion
        }

        // KERNELS
        static void GPU_IntegrationKernal(Index1 index, ArrayView<float> OutPut, float dz, float Ogamm0, float Om0, float Ode0, float neff_per_nu, float nmasslessnu, float nu_y)
        {
            float opz = 1f + (dz * index);
            float k = 0.3173f / opz;
            float Or0 = (Ogamm0 * (1f + 0.22710731766f * neff_per_nu * (nmasslessnu + XMath.Pow(1f + XMath.Pow(k * nu_y, 1.83f), 0.54644808743f))));

            OutPut[index] = XMath.Rsqrt(XMath.Pow(opz, 3f) * (opz * Or0 + Om0) + Ode0) * dz; 
        }


    }

}
