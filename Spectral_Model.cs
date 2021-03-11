using ILGPU;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MachineLearningSpectralFittingCode
{
    class Spectral_Model
    {
        public Spectral_Model(string path,
                       bool milky_Way_Reddening = true,
                       bool hPF_Mode = true,
                       ushort n_Masked_Amstrongs = 20)
        {
            this.Path = path;
            this.Milky_Way_Reddening = milky_Way_Reddening;
            this.HPF_Mode = hPF_Mode;
            this.N_Masked_Amstrongs = n_Masked_Amstrongs;
        }

        // SPECTRA SPECIFIC VARIABLES
        #region
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
        public Vector Instrument_Resolution { get; private set; }

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
        #endregion

        // MODEL SPECIFIC VARIABLES
        #region
        int velocity_dispersion_r { get; set; }
        int fit_per_iteration_cap { get; set; }
        List<double> delta_lamdba_lib { get; set; }
        double deltal { get; set; }


        #endregion

        // Spectrum Initialisation
        public async void InitialiseSpectraParameters(Accelerator gpu, Vector Data, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
        {
            bool warn = false;

        retryWave:
            try
            {
                // Set Long tasks to execute
                this.Wavelength = await Vector.AccessSliceAsync(gpu, Data, 0, 'c');  // WAVELENGTH
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryWave;
            }

        retryFlux:
            try
            {
                // Set Long tasks to execute
                this.Flux = await Vector.AccessSliceAsync(gpu, Data, 1, 'c');  // FLUX
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryFlux;
            }

        retryErr:
            try
            {
                // Set Long tasks to execute
                this.Error = await Vector.AccessSliceAsync(gpu, Data, 2, 'c');  // ERROR
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryErr;
            }


            this.Redshift = Redshift;
            this.RA_DEC = RA_DEC;
            this.Velocity_Dispersion = Velocity_Disp;
            this.Trust_Flag = 1;
            this.ObjID = 0;

        retryInstru:
            try
            {
                this.Instrument_Resolution = Vector.Fill(gpu, instrument_resolution, this.Wavelength.Value.Length);
            }
            catch (Exception)
            {
                warn = true;
                Task.Delay(100);
                goto retryInstru;
            }


            if (this.Milky_Way_Reddening)
            {
                this.ebv_MW = this.GetDustRADEC(this.RA_DEC, "ebv");
            }
            else
            {
                this.ebv_MW = 0f;
            }


        retryRestWave:
            try
            {
                this.Restframe_Wavelength = Vector.ScalarOperation(gpu, this.Wavelength, (1 + this.Redshift), '/');
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryRestWave;
            }

        // Remove Bad data from the spectrum
        // Generates the Filter Mask
        retryBadDat:
            try
            {
                float[] BadDataMask = this.GenerateDataMask(gpu, this.Flux, this.Error);

                if (BadDataMask.Contains(0f))
                {
                    // Filter Out the bad data
                    Console.WriteLine("Warning Bad Data Detected");
                }
                else
                {
                    //Console.WriteLine("Data is Fine");
                }
                // Else just proceed as the Data is Fine   
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryBadDat;
            }



        retryDistLum:
            try
            {
                // CALCULATE LUMINOSITY DISTANCE in CM
                // Two versions to Test Normal and Optimised Algo. no noticable performance difference noted so far
                this.Distance_Luminosity = UtilityMethods.Mpc2cm(await Program.cosmology.GPU_IntegrationOptiAsync(gpu, this.Redshift, 1e-8f)); //Program.cosmology.luminosity_distance(this.Redshift));
            }
            catch (Exception)
            {
                warn = true;
                await Task.Delay(100);
                goto retryDistLum;
            }

            
            // Initialise Model Section
            this.InitialiseSPModel();


            // Display warning
            if (warn)
            {
                Console.WriteLine("Warning Code Executed Slower Than Expected - Insufficient Memory");
            }

        }

        // Model Initialisation
        private void InitialiseSPModel()
        {
            this.delta_lamdba_lib = new List<double>();
            switch (Program.config.Model_Key)
            {
                case 0b00001_001:
                    this.delta_lamdba_lib.Add(2.55d);
                    break;
                case 0b00010_001:
                    this.delta_lamdba_lib.Add(3.40d);
                    break;
                case 0b00100_001:
                    this.delta_lamdba_lib.Add(0.55d);
                    break;
                case 0b01000_001:
                    this.delta_lamdba_lib.Add(0.10d);
                    break;


                case 0b00001_010:
                    this.delta_lamdba_lib.AddRange(Constants.r_model);
                    break;
                case 0b00010_010:
                    this.delta_lamdba_lib.AddRange(Constants.r_model);
                    break;

                default:
                    throw new Exception("Incorrect Model key");
            }
            this.velocity_dispersion_r = (int)(MathF.Round(this.Velocity_Dispersion / 5f) * 5f);
        }


        // SPECTRAL - SPECTRUM FUNCTIONS
        #region
        private float[] GenerateDataMask(Accelerator gpu, Vector flux, Vector Error) // Also add emission lines filter
        {

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(GPU_GenerateDataMaskKernal);

            var buffer = gpu.Allocate<float>(flux.Value.Length);
            var buffer2 = gpu.Allocate<float>(flux.Value.Length);
            var buffer3 = gpu.Allocate<float>(Error.Value.Length);

            buffer.MemSetToZero();
            buffer2.MemSetToZero();
            buffer3.MemSetToZero();

            buffer2.CopyFrom(flux.Value, 0, 0, flux.Value.Length);
            buffer3.CopyFrom(Error.Value, 0, 0, Error.Value.Length);


            kernelWithStream(Stream, buffer.Length, buffer.View, buffer2.View, buffer3.View);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

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

        #endregion

        // MODELS - MODEL FUNCTIONS
        #region
        public void Fit_models_to_data()
        {
            // enumerate over each model_lib (mi, mm)?
            // only 1 model library so, index mi = 0, mm = "MaStar"


            // model_lib = Model_key = MaStar-Th
            // deltal_libs = List<double> of r_model


            // enumerate over imfs (ii)
            // only 1 imf so, foreach ii in imfs
            this.deltal = this.delta_lamdba_lib[0];
            this.Get_Model(); //sets: model_wave_int,model_flux_int,age,metal


            // this.raw_model_wave_int = model_wave_int
            // this.raw_model_flux_int ... this.metal = metal
            Match_data_models();
            normalise_spec();

            // Part 3 - correction from dust attenuation
            // -- LOG TIME COMPARISON -- CHECKPOINT
            if (Program.config.HPF_Mode)
            {

            }
            else // hpf only if False
            {

            }


            // Part 5
            // Gets the best model

        }

        private void Get_Model()
        {
            if (Program.config.Model_Key % 2 == 1)
            {
                //var first_file = true;
                // List of model_files = []

                // Constructs the metallicity array of models:
                //var all_metal_files = sorted(glob.glob(model_path+"*"))
                // metal_files = []
                // metal = [] #[-2.25,-1.35,-0.33,0,0.35]
                // for z in range(len(all_metal_files)):
                //.
                //.
                //.


            }
            else if (Program.config.Model_Key % 2 == 0)
            {
                // This section can go to the config pre initialisation
                var slope = 0f;
                switch (Program.config.IMF)
                {
                    case 0: // kr
                        slope = 1.3f;
                        break;
                    case 1: // ss
                        slope = 2.35f;
                        break;

                    default:
                        throw new Exception("Unrecognised IMF");
                }
                var sidx = Array.IndexOf(Constants.s, slope);
                // This section - END

                List<float> model_flux = new List<float>();
                List<float> age_model = new List<float>();
                List<float> metal_model = new List<float>();

                for (int i = 0; i < Constants.t.Length; i++)
                {
                    //index i in t , t value
                    if (Constants.t[i] < Program.config.MinMax_Model_Age[0] || Constants.t[i] > Program.config.MinMax_Model_Age[1])
                    {
                        continue;
                    }

                    for (int j = 0; j < Constants.Z.Length; j++)
                    {
                        // index, Z value
                        if (Constants.Z[j] < Program.config.MinMax_Model_Metalicity[0] || Constants.Z[j] > Program.config.MinMax_Model_Metalicity[1])
                        {
                            continue;
                        }
                        if (Constants.Z[j] < -1.35f && Constants.t[i] < 1)
                        {
                            continue;
                        }

                        float[] flux = new float[4563];
                        for (int k = 0; k < 4563; k++)
                        {
                            flux[k] = Constants.fluxgrid[i, j, sidx, k];
                        }

                        // no conversion to vacuum needed, assuming models are in vacuum
                        
                        // downgrades the model
                        if (Program.config.Downgrade_models)
                        {
                            flux = DownGrade(Constants.wavelength, flux, this.deltal, this.velocity_dispersion_r,this.Restframe_Wavelength, this.Instrument_Resolution);
                        }
                        else
                        {
                            // flux = flux
                        }

                        // Reddens the models
                        if (this.ebv_MW != 0)
                        {
                            var attenuations = unred(Constants.wavelength, -this.ebv_MW); // ebv = 0f - ebv_mw

                            try
                            {
                                model_flux.AddRange(Vector.ConsecutiveProduct(Program.gpu, (new Vector(flux)), (new Vector(attenuations))).Value);

                            }
                            catch (Exception)
                            {
                                throw new Exception("Spectral Model Run Out of Memory In Fit_Data_To_Models - Get_Models Function");
                            }

                        }
                        else
                        {
                            model_flux.AddRange(flux);
                        }

                        age_model.Add(i);
                        metal_model.Add(MathF.Pow(10f, j));

                        
                    }

                }

            }
            /*
             * sets the
             * model_wave_int
             * model_flux_int
             * age
             * metal
             */

        }

        private void Match_data_models()
        {
            /*
             * sets the
             * wave
             * model_flux_raw
             * 
             * outputs the 
             * data_flux
             * error_flux
             */
        }

        private void normalise_spec()
        {
            /*
             * sets the 
             * model_flux
             * mass_factors
             */
        }

        private float[] DownGrade(float[] wavelength, float[] flux, double deltal, int vdisp_round, Vector rest_wavelength, Vector r_instrument)
        {


            return flux;
        }

        private float[] unred(float[] wavelength, float ebv_mw)
        {


            return wavelength;
        }


        #endregion

    }



}
