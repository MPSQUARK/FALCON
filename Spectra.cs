using ILGPU;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    public class Spectra
    {
        public Spectra(string path,
                       bool milky_Way_Reddening = true,
                       bool hPF_Mode = true,
                       ushort n_Masked_Amstrongs = 20)
        {
            this.Path = path;
            this.Milky_Way_Reddening = milky_Way_Reddening;
            this.HPF_Mode = hPF_Mode;
            this.N_Masked_Amstrongs = n_Masked_Amstrongs;
        }

        // constructor variables
        public string Path { get; private set; }
        public bool Milky_Way_Reddening { get; private set; }
        public bool HPF_Mode { get; private set; }
        public ushort N_Masked_Amstrongs { get; private set; }

        // Initialisation variables
        public Vector Wavelength { get; private set; }
        public Vector Flux { get; private set; }
        public Vector Error { get; private set; }
        public float Redshift { get; private set; }
        public float[] RA_DEC { get; private set; }
        public float Velocity_Dispersion { get; private set; }
        public float Instrument_Resolution { get; private set; }

        public float Distance_Luminosity { get; private set; }
        public Vector Bad_Flags { get; private set; }
        public Vector Restframe_Wavelength { get; private set; }
        // !!errorNES!!
        public byte Trust_Flag { get; private set; }
        // !!errorNES!!
        public byte ObjID { get; private set; }
        // Bad Data ?? 
        // mask emission lines
        public float ebv_MW { get; private set; }




        public void InitialiseSpectraParameters(AcceleratorId acceleratorId, Vector Data, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
        {
            this.Wavelength = Vector.AccessSlice(acceleratorId, Data, 0, 'c');  // GPU Batch 1
            this.Flux = Vector.AccessSlice(acceleratorId, Data, 1, 'c');        // GPU Batch 1
            this.Error = Vector.AccessSlice(acceleratorId, Data, 2, 'c');       // GPU Batch 1
            this.Redshift = Redshift;
            this.RA_DEC = RA_DEC;
            this.Velocity_Dispersion = Velocity_Disp;
            this.Instrument_Resolution = instrument_resolution;

            // CALCULATE LUMINOSITY DISTANCE in CM
            this.Distance_Luminosity = UtilityMethods.Mpc2cm(Program.cosmology.luminosity_distance(this.Redshift)); // GPU Batch 1

            //this.Bad_Flags = Vector.Fill(acceleratorId, 1, this.Wavelength.Value.Length); // Unnessessary 
            this.Restframe_Wavelength = Vector.ScalarOperation(acceleratorId, this.Wavelength, (1 + this.Redshift), '/'); // GPU Batch 2
            this.Trust_Flag = 1;
            this.ObjID = 0;


            // Remove Bad data from the spectrum

            // Generates the Filter Mask
            float[] BadDataMask = this.GenerateDataMask(acceleratorId, this.Flux, this.Error); // GPU Batch 2


            if (BadDataMask.Contains(0f))
            {
                // Filter Out the bad data
                Console.WriteLine("Warning Bad Data Detected");
            }
            else
            {
                Console.WriteLine("Data is Fine");
            }

            // Else just proceed as the Data is Fine

            if (this.Milky_Way_Reddening)
            {
                this.ebv_MW = this.GetDustRADEC(this.RA_DEC, "ebv");
            }
            else
            {
                this.ebv_MW = 0f;
            }

        }

        private float[] GenerateDataMask(AcceleratorId acceleratorId, Vector flux, Vector Error) // Also add emission lines filter
        {

            using var context = new Context();
            context.EnableAlgorithms();

            using var accelerator = Accelerator.Create(context, acceleratorId);

            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(GPU_GenerateDataMaskKernal);

            var buffer = accelerator.Allocate<float>(flux.Value.Length);
            var buffer2 = accelerator.Allocate<float>(flux.Value.Length);
            var buffer3 = accelerator.Allocate<float>(Error.Value.Length);
            
            buffer.MemSetToZero();
            buffer2.MemSetToZero();
            buffer3.MemSetToZero();

            buffer2.CopyFrom(flux.Value, 0, 0, flux.Value.Length);
            buffer3.CopyFrom(Error.Value, 0, 0, Error.Value.Length);


            kernel(buffer.Length, buffer.View, buffer2.View, buffer3.View);

            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();

            return Output;
        }

        // GPU KERNEL
        static void GPU_GenerateDataMaskKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> flux, ArrayView<float> error)
        {
            // False means Exclude/Bad Data, True means Good Data
            OutPut[index] = Convert.ToSingle(!(float.IsNaN(flux[index]) || float.IsInfinity(flux[index]) || (flux[index] <= 0f) || float.IsNaN(error[index]) || float.IsInfinity(error[index])));
        }




        private float GetDustRADEC(float[] RADEC, string dustmap, bool interpolate = true)
        {

            float[] l_b = Eq2Gal(
                UtilityMethods.Degree2Radians(RADEC[0]), 
                UtilityMethods.Degree2Radians(RADEC[1])
                );

            return Get_SFD_dust(l_b[0], l_b[1]);
        }

        private float[] Eq2Gal(float ra, float dec)
        {

            float Cos_dec_Mult_Cos_GNEV0 = (float)(Math.Cos(dec) * Math.Cos(ra - Constants.Galactic_Northpole_Equatorial.Value[0]));

            float b = (float)Math.Asin(
                Math.Sin(dec) * Math.Sin(Constants.Galactic_Northpole_Equatorial.Value[1]) + 
                Math.Cos(Constants.Galactic_Northpole_Equatorial.Value[1]) * Cos_dec_Mult_Cos_GNEV0
                );

            float l = (float)Math.Atan2(
                Math.Sin(dec) * Math.Cos(Constants.Galactic_Northpole_Equatorial.Value[1]) - 
                Math.Sin(Constants.Galactic_Northpole_Equatorial.Value[1]) * Cos_dec_Mult_Cos_GNEV0, 
                Math.Cos(dec) * Math.Sin(ra - Constants.Galactic_Northpole_Equatorial.Value[0])
                ) + UtilityMethods.Degree2Radians(33f);

            if (l < 0)
            {
                l += Constants.TwoPi;
            }

            l %= Constants.TwoPi;


            return new float[2] { l * Constants.Rad2DegFactor, b * Constants.Rad2DegFactor };
        }

        private float Get_SFD_dust(float logitude, float latitude, string dustmap = "ebv", bool interpolate = true)
        {
            // get path to dust map



            return 0f;
        }

    }
}
