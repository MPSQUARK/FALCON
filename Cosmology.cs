using ILGPU;
using ILGPU.Runtime;
using ILGPU.Algorithms;
using BAVCL;

namespace FALCON;

public class Cosmology
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
    readonly bool flat = true;
    readonly float[] m_nu;
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

    private Func<float, Inv_efunc_scalar_args_struct, float> inv_efunc_scalar;
    private Inv_efunc_scalar_args_struct inv_efunc_scalar_args;

    // optimisation parameters
    float
        Or0_OptiA,
        Or0_OptiB,
        nuyp_Opti;


    public Cosmology()
    {
        Oc0 = 0.2589f;
        Ob0 = 0.04860f;
        Om0 = 0.3075f;
        H0 = 67.74f;
        n = 0.9667f;
        sigma8 = 0.8159f;
        tau = 0.066f;
        z_reion = 8.8f;
        t0 = 13.799f;
        Tcmb0 = 2.7255f;
        Neff = 3.046f;
        flat = true;
        m_nu = [0f, 0f, 0.06f];

        Initialise();
    }

    public void Initialise()
    {
        Odm0 = Om0 - Ob0;
        h = H0 * 0.01f;
        hubble_distance = (Constants.c * 0.001f) / H0;
        H0_s = H0 * Constants.H0units_to_invs;
        hubble_time = (float)(Constants.sec_to_Gyr / H0_s);
        critical_density0 = Constants.critdens_const * Math.Pow(H0_s, 2);
        nneutrinos = (int)MathF.Floor(Neff);

        massivenu = false;
        if (nneutrinos > 0 && Tcmb0 > 0)
        {
            neff_per_nu = Neff / nneutrinos;

            if (m_nu.Max() == 0)
            {
                nmasslessnu = nneutrinos;
                nmassivenu = 0;
            }
            else
            {
                massivenu = true;
                if (m_nu.Length != nneutrinos)
                {
                    throw new Exception("Class Cosmology : func Initialise - Unexpected number of neutrino masses");
                }

                massivenu_mass = (from nu in m_nu
                                  where nu != 0
                                  select nu).ToArray();
                nmassivenu = massivenu_mass.Length;
                nmasslessnu = nneutrinos - nmassivenu;
            }
        }

        if (Tcmb0 > 0)
        {
            Ogamma0 = (Constants.a_B_c2 * Math.Pow(Tcmb0, 4d) / critical_density0);
            Tnu0 = 0.7137658555036082f * Tcmb0;

            if (massivenu)
            {
                int length = massivenu_mass.Length;
                nu_y = new float[length];

                for (int i = 0; i < length; i++)
                {
                    nu_y[i] = (float)(massivenu_mass[i] / (Constants.kB_evK * Tnu0));
                }


                Onu0 = (float)(Ogamma0 * Nu_relative_density(0f));

            }
            else
            {
                Onu0 = (float)(0.22710731766f * Neff * Ogamma0);
            }
        }
        else
        {
            throw new Exception("Class Cosmology : func Initialise - Unexpected value of Tcmb0, should be tcmb0 > 0 Kelvin");
        }

        Ode0 = (float)(1f - (Om0 + Ogamma0 + Onu0));
        Ok0 = 0; // Assuming no curvature

        if (!massivenu)
        {
            inv_efunc_scalar = Flcdm_inv_efunc_nomnu;
            inv_efunc_scalar_args = new Inv_efunc_scalar_args_struct
            {
                Om0 = Om0,
                Ode0 = Ode0,
                Or0 = (float)(Ogamma0 + Onu0)
            };
        }
        else
        {
            inv_efunc_scalar = Flcdm_inv_efunc;
            inv_efunc_scalar_args = new Inv_efunc_scalar_args_struct
            {
                Om0 = Om0,
                Ode0 = Ode0,
                Ogamma0 = (float)Ogamma0,
                neff_per_nu = neff_per_nu,
                nmasslessnu = nmasslessnu,
                nu_y = nu_y
            };
        }


        // Further optimisations
        Or0_OptiA = (float)(Ogamma0 + (0.22710731766f * neff_per_nu * Ogamma0 * nmasslessnu));
        Or0_OptiB = (float)(0.22710731766f * neff_per_nu * Ogamma0);
        nuyp_Opti = MathF.Pow(nu_y[0], 1.83f);

    }


    private struct Inv_efunc_scalar_args_struct
    {
        public float Om0;
        public float Ode0;
        public float Or0;

        public float Ogamma0;
        public float neff_per_nu;
        public float nmasslessnu;
        public float[] nu_y;
    }


    private float Nu_relative_density(float z)
    {
        // See Komatsu et al. 2011, eq 26 and the surrounding discussion
        float prefac = 0.22710731766f;  // 7/8 (4/11)^4/3 -- see any cosmo book

        // Assume that z is always a scalar
        if (!massivenu)
        {
            return prefac * Neff;
        }

        // These are purely fitting constants -- see the Komatsu paper
        float p = 1.83f;
        float invp = 0.54644808743f;  // 1.0 / p
        float k = 0.3173f;

        float[] rel_mass_per = new float[nu_y.Length];
        for (int i = 0; i < nu_y.Length; i++)
        {
            rel_mass_per[i] = MathF.Pow(1f + MathF.Pow(k * (nu_y[i] / (1f + z)), p), invp);
        }
        float rel_mass = rel_mass_per.Sum() + nmasslessnu;

        return prefac * neff_per_nu * rel_mass;
    }

    private float Flcdm_inv_efunc_nomnu(float z, Inv_efunc_scalar_args_struct args)
    {
        float opz = 1 + z;
        return MathF.Pow(MathF.Pow(opz, 3f) * (opz * args.Or0 + args.Om0) + args.Ode0, -0.5f);
    }

    private float Flcdm_inv_efunc(float z, Inv_efunc_scalar_args_struct args)
    {
        float opz = 1f + z;
        float Or0 = args.Ogamma0 * (1f + Nufunc(opz, args.neff_per_nu, args.nmasslessnu, args.nu_y));
        return MathF.Pow(MathF.Pow(opz, 3f) * (opz * Or0 + args.Om0) + args.Ode0, -0.5f);
    }

    private static float Nufunc(float opz, float neff_per_nu, float nmasslessnu, float[] nu_y)
    {
        int N = nu_y.Length;
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

    public float Luminosity_distance(GPU gpu, float redshift)
    {
        return (1f + redshift) * Comoving_transverse_distance(gpu, redshift);
    }

    public float Comoving_transverse_distance(GPU gpu, float redshift)
    {
        float dc = Integral_comoving_distance(gpu, redshift);
        return dc;
    }

    public float Integral_comoving_distance(GPU gpu, float redshift)
    {

        return hubble_distance * GPU_Integration(gpu, redshift, 1e-8f);

    }


    private static float Integrate(Func<float, Inv_efunc_scalar_args_struct, float> func, float z, float da, Inv_efunc_scalar_args_struct args)
    {
        int itter = (int)(z / da);
        float[] vals = new float[itter];

        for (int i = 0; i < itter; i++)
        {
            vals[i] = func(i * da, args) * da;
        }

        return vals.Sum();
    }


    //// Integrates the cosmology func to get luminosity distance 
    //public float GPU_Integration(Accelerator gpu, float z, float dz)
    //{
    //    int length = (int)(z / dz);

    //    AcceleratorStream Stream = gpu.CreateStream();

    //    var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, float, float, float, float, float, float, float>(GPU_IntegrationKernal);

    //    var buffer = gpu.Allocate<float>(length);
    //    buffer.MemSetToZero(Stream);

    //    kernelWithStream(Stream, buffer.Length, buffer.View, dz, (float)this.Ogamma0, this.Om0, this.Ode0, this.neff_per_nu, this.nmasslessnu, this.nu_y[0]);

    //    Stream.Synchronize();

    //    float[] Output = buffer.GetAsArray(Stream);

    //    buffer.Dispose();

    //    Stream.Dispose();

    //    return (1f + z) * this.hubble_distance * Output.Sum();
    //}


    // Integrates the cosmology func to get luminosity distance Asynchronously
    public float GPU_Integration(GPU gpu, float z, float dz)
    {
        // WARNING THIS NEEDS ATTENTION!!!
        Vector opz = Vector.Linspace(gpu, 1, z, XMath.Abs((int)(z / dz)));
        // WARNING THIS NEEDS ATTENTION!!!

        opz.IncrementLiveCount();
        MemoryBuffer1D<float, Stride1D.Dense> opzBuffer = opz.GetBuffer();

        var kernel = gpu.accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, float, float, float, float, float, float, float>(GPU_IntegrationKernal);
        kernel(gpu.accelerator.DefaultStream, opzBuffer.IntExtent, opzBuffer.View, dz, (float)Ogamma0, Om0, Ode0, neff_per_nu, nmasslessnu, nu_y[0]);
        gpu.accelerator.Synchronize();

        float sum = opz.Sum();

        opz.DecrementLiveCount();

        // WARNING THIS NEEDS ATTENTION!!!
        if (z < 0) { sum = -sum; }
        // WARNING THIS NEEDS ATTENTION!!!


        return (1f + z) * hubble_distance * sum;
    }


    // KERNELS
    static void GPU_IntegrationKernal(Index1D index, ArrayView<float> OutPut, float dz, float Ogamma0, float Om0, float Ode0, float neff_per_nu, float nmasslessnu, float nu_y)
    {
        float opz = 1f + (dz * index);
        float k = 0.3173f / opz;
        float Or0 = (Ogamma0 * (1f + 0.22710731766f * neff_per_nu * (nmasslessnu + XMath.Pow(1f + XMath.Pow(k * nu_y, 1.83f), 0.54644808743f))));

        OutPut[index] = XMath.Rsqrt(XMath.Pow(opz, 3f) * (opz * Or0 + Om0) + Ode0) * dz;
    }

}