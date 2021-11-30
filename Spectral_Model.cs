using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BAVCL;
using System.Runtime.CompilerServices;

namespace FALCON
{
    class Spectral_Model
    {
        public Spectral_Model(string path,
                       bool milky_Way_Reddening = true, // SHOULD BE TRUE
                       bool hPF_Mode = true,
                       ushort n_Masked_Amstrongs = 20,
                       GPU GPU = null)
        {
            this.Path = path;
            this.Milky_Way_Reddening = milky_Way_Reddening;
            this.HPF_Mode = hPF_Mode;
            this.N_Masked_Amstrongs = n_Masked_Amstrongs;
            this.gpu = GPU;
        }

        GPU gpu;

        // SPECTRA SPECIFIC VARIABLES
        #region
        // constructor variables
        public string Path { get; private set; }
        public bool Milky_Way_Reddening { get; private set; }
        public bool HPF_Mode { get; private set; }
        public ushort N_Masked_Amstrongs { get; private set; }

        // Initialisation variables
        /// <summary>
        /// DATA wavelength in OBSERVERS Frame of REFERENCE
        /// </summary>
        public Vector Wavelength { get; private set; }
        /// <summary>
        /// DATA Flux Parameter
        /// </summary>
        public Vector Flux { get; private set; }
        /// <summary>
        /// DATA ERROR PARAMETER
        /// </summary>
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
        public bool DataLiesOutsideModel = false;
        #endregion

        // MODEL SPECIFIC VARIABLES
        #region
        int velocity_dispersion_r { get; set; }
        int fit_per_iteration_cap { get; set; }
        // Model Main Values
        /// <summary>
        /// MODEL wavelengths
        /// </summary>
        public Vector Model_wavelength;
        /// <summary>
        /// MODEL FLUXES - List of MODEL FLUXES
        /// </summary>
        public Vector Model_flux; // Becomes NORMED after Method - NormaliseSpec
        /// <summary>
        /// MODEL AGES
        /// </summary>
        public float[] Model_ages;
        /// <summary>
        /// MODEL METALICITIES
        /// </summary>
        public float[] Model_metals;
        public float[] Mass_factor;
        public int[] MatchingWavelengthMapping;


        #endregion

        // Spectrum Initialisation
        public void InitialiseSpectraParameters(BAVCL.Vector Data, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
        {

            this.Wavelength = Data.GetColumnAsVector(0);
            this.Flux = Data.GetColumnAsVector(1);
            this.Error = Data.GetColumnAsVector(2);

            this.Redshift = Redshift;
            this.RA_DEC = RA_DEC;
            this.Velocity_Dispersion = Velocity_Disp;
            this.Trust_Flag = 1;
            this.ObjID = 0;

            this.Instrument_Resolution = Vector.Fill(gpu, instrument_resolution, this.Wavelength.Value.Length);

            if (this.Milky_Way_Reddening)
            {
                this.ebv_MW = this.GetDustRADEC(this.RA_DEC, "ebv");
            }
            else
            {
                this.ebv_MW = 0f;
            }

            this.Restframe_Wavelength = this.Wavelength * (1 / (1 + this.Redshift));

            // Generates the Filter Mask 2 Remove Bad data from the spectrum
            Vector BadDataMask = this.GenerateDataMask(this.Flux, this.Error);
            BadDataMask.SyncCPU();

            if (BadDataMask.Value.Contains(0f))
            {
                // Filter Out the bad data
                Console.WriteLine("Warning Bad Data Detected");
            }
            else
            {
                //Console.WriteLine("Data is Fine");
            } 

            // CALCULATE LUMINOSITY DISTANCE in CM and convert to Mpc
            this.Distance_Luminosity = UtilityMethods.Mpc2cm(Program.cosmology.GPU_Integration(gpu, this.Redshift, 1e-8f));


            // Initialise Model Section
            this.InitialiseSPModel();

        }

        public void InitialiseSpectraParameters(Vector Wavelength, Vector Flux, Vector Error, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
        {

            // Set Long tasks to execute
            this.Wavelength = Wavelength;

            // Set Long tasks to execute
            this.Flux = Flux;

            // Set Long tasks to execute
            this.Error = Error;

            this.Redshift = Redshift;
            this.RA_DEC = RA_DEC;
            this.Velocity_Dispersion = Velocity_Disp;
            this.Trust_Flag = 1;
            this.ObjID = 0;


            if (this.Milky_Way_Reddening)
            {
                this.ebv_MW = this.GetDustRADEC(this.RA_DEC, "ebv");
            }
            else
            {
                this.ebv_MW = 0f;
            }

            this.Restframe_Wavelength = this.Wavelength * (1f / (1+this.Redshift));


            // Remove Bad data from the spectrum
            // Generates the Filter Mask

            Vector BadDataMask = this.GenerateDataMask(this.Flux, this.Error);
            BadDataMask.SyncCPU();

            if (BadDataMask.Value.Contains(0f))
            {
                // Filter Out the bad data
                int goodvals = (int)BadDataMask.Sum();
                    
                float[] new_wavelength = new float[goodvals];
                float[] new_restframewavelength = new float[goodvals];
                float[] new_flux = new float[goodvals];
                float[] new_error = new float[goodvals];

                for (int i = 0, j = 0; i < BadDataMask.Length; i++)
                {
                    if (BadDataMask.Value[i] != 1) { continue; }

                    new_wavelength[j] = this.Wavelength.Value[i];
                    new_restframewavelength[j] = this.Restframe_Wavelength.Value[i];
                    new_flux[j] = this.Flux.Value[i];
                    new_error[j] = this.Error.Value[i];
                    j++;
                }

                this.Wavelength.UpdateCache(new_wavelength);
                this.Restframe_Wavelength.UpdateCache(new_restframewavelength);
                this.Flux.UpdateCache(new_flux);
                this.Error.UpdateCache(new_error);

                //Console.WriteLine($"Warning Bad Data Detected, {BadDataMask.Sum() / BadDataMask.Length * 100}% are Good");
            }
            else
            {
                //Console.WriteLine("Data is Fine");
            }

            this.Instrument_Resolution = Vector.Fill(gpu, instrument_resolution, this.Wavelength.Length);

            // CALCULATE LUMINOSITY DISTANCE in CM and convert to Mpc
            this.Distance_Luminosity = UtilityMethods.Mpc2cm(Program.cosmology.GPU_Integration(gpu, this.Redshift, 1e-8f));


            // Initialise Model Section
            this.InitialiseSPModel();

        }



        // Model Initialisation
        private void InitialiseSPModel()
        {
            this.velocity_dispersion_r = (int)(MathF.Round(this.Velocity_Dispersion / 5f) * 5f);
        }


        // SPECTRAL - SPECTRUM FUNCTIONS
        #region
        private Vector GenerateDataMask(Vector flux, Vector error) // Also add emission lines filter
        {
            flux.IncrementLiveCount();
            error.IncrementLiveCount();

            Vector output = new Vector(gpu, new float[flux.Value.Length]);
            output.IncrementLiveCount();

            MemoryBuffer<float> buffer = output.GetBuffer();
            MemoryBuffer<float> buffer2 = flux.GetBuffer();
            MemoryBuffer<float> buffer3 = error.GetBuffer();

            var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(GPU_GenerateDataMaskKernal);
            kernelWithStream(gpu.accelerator.DefaultStream, buffer.Length, buffer.View, buffer2.View, buffer3.View);

            gpu.accelerator.Synchronize();

            flux.DecrementLiveCount();
            Error.DecrementLiveCount();
            output.DecrementLiveCount();

            return output;
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

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private float[] Eq2Gal(float ra, float dec)
        {

            float Cos_dec_Mult_Cos_GNEV0 = (float)(Math.Cos(dec) * Math.Cos(ra - Constants.Galactic_Northpole_Equatorial[0]));

            float b = (float)Math.Asin(
                Math.Sin(dec) * Math.Sin(Constants.Galactic_Northpole_Equatorial[1]) +
                Math.Cos(Constants.Galactic_Northpole_Equatorial[1]) * Cos_dec_Mult_Cos_GNEV0
                );

            float l = (float)Math.Atan2(
                Math.Sin(dec) * Math.Cos(Constants.Galactic_Northpole_Equatorial[1]) -
                Math.Sin(Constants.Galactic_Northpole_Equatorial[1]) * Cos_dec_Mult_Cos_GNEV0,
                Math.Cos(dec) * Math.Sin(ra - Constants.Galactic_Northpole_Equatorial[0])
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
            //Console.WriteLine($"longitude {logitude} , latitude {latitude}");

            // get path to dust map

            // if var dustmap is NOT a string 
            //      - Raise Exception
            //  + ALWAYS WILL BE A STRING

            // var dml = dustmap in lowercase

            /* if dustmap is 'ebv' 'eb-b' or 'e(b-v)'
             *  + ONLY accepted variation will be 'ebv'
             * 
             *      - read in file at os.environ['FF_DIR']+'/data/SFD_dust_4096_%s.fits'
             * + TEST WHERE THIS IS TO CONVERT TO HDF5
             * 
             * + RETURN ALL OTHER TYPES AS DUSTMAP NOT FOUND
             * 
             * if dustmap = i100
             * if dustmap = x
             * if dustmap = t
             * if dustmap = mask
             * else dustmapfn = dustmap ?? would that not crash?
             * 
             * if longitude is a scalar
             *      - l = convert longitude to an array and multipy by pi/180
             * + WHICH it IS
             * if latitude is a scalar
             *      - b = convert latitude to an array and multipy by pi/180
             *      
             * if length of l and b DONT match, RAISE EXCEPTION
             *  + WILL ALWAYS MATCH SINCE l and b are floats
             * 
             * if NOT dustmapfn.contains(%s)
             *  + %s WILL ALWAYS be in the path name
             *      - f = fits.open(dustmapfn)
             *      - try : mapds = [f[0].data] then f.close()
             * + READ IN DUSTMAP in python and convert to HDF5 for this and read in PreInitialise Phase
             * 
             *      - assert that mapds[-1].shape[0] is = to mapds[-1].shape[1] 
             * + Test in Python what does means
             * 
             *      - polename = dustmapfn.split('.')[0].split('_')[-1].lower() ?? Can do a Contains or Regex test instead??
             * 
             *      - if polename is 'ngp'
             *          -> n=[1]
             *      - if sum( b > 0 ) > 0
             *          -> b = b
             *      - if polename is 'sgp'
             *          -> n = [-1]
             *          - if sum ( b < 0 ) > 0:
             *              -> b = b
             *      - else
             *          - RAISE EXCEPTION
             *      
             *      -> masks = [ones(same shape as b) converted to boolean ]
             * ELSE
             *      -> nmask = b >= 0 
             *      -> smask = !nmask
             *      
             *      masks = [nmask, smask]
             *      ns = [1,-1]
             *      
             *      mapds = []
             *     =>-> READ IN THE dustmapfn fits for 'ngp'/'sgp'
             *      try: mapds.add(f[0].data) then close(f)
             *      asstert that mapds[-1].shape[0] is EQUAL to mapds[-1].shape[1]
             *     <=->repeat for 'sgp'
             *      
             *      retvals = []
             *      for n,mapd,m in zip(ns,mapds,masks) 
             *          + CHECK THIS FUNCTIONALITY IN PYTHON TO REPLICATE
             *          npix = mapd.shape[0]
             *      
             *      
             *      
            */



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
            this.Get_Model(); //sets: model_wave_int,model_flux_int,age,metal


            // this.raw_model_wave_int = model_wave_int
            // this.raw_model_flux_int ... this.metal = metal
            //Match_data_models(); 
            NormaliseSpec();

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


            // TEST CODE
            GenerateInterpSSP();
            


        }

        private void Get_Model()
        {

            // M11 Models
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

                return;
            }
            
            // MaStar Models
            if (Program.config.Model_Key % 2 == 0)
            {
                // Gets Wavelengths from Model Data
                this.Model_wavelength = Constants.wavelength;  // Lacking 1dp precision
                // Gets Indexes of Models matching wavlengths
                this.TrimModel();

                (int[] indx, Spectral_resolution res, VariableGaussianKernel vGK) = this.DownGradeModelInvar(this.Model_wavelength, this.velocity_dispersion_r, Constants.sres, this.Instrument_Resolution);

                List<float> model_flux = new List<float>();
                List<float> age_model = new List<float>();
                List<float> metal_model = new List<float>();


                Vector flux = new Vector(gpu, new float[Constants.wavelength.Length], 1, false);
                

                // MAX AGE needs fixing to use data from the fits input

                for (int i = 0; i < Constants.t.Length; i++)
                {
                    //index i in t , t value
                    if ((Constants.t[i] < Program.config.MinMax_Model_Age[0]) || (Constants.t[i] > Program.config.MinMax_Model_Age[1]))
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

                        // TRIM DOWN THE WAVELENGTH AND THE FLUX to correspond to Data
                        for (int k = 0; k < flux.Length; k++)
                        {
                            flux.Value[k] = Constants.fluxgrid[i, j, Constants.sidx, k];
                        }
                        flux.UpdateCache();

                        // no conversion to vacuum needed, assuming models are in vacuum

                        // downgrades the model
                        if (Program.config.Downgrade_models)
                        {
                            flux = this.DownGrade(gpu,flux, indx, res, vGK);
                        }
                        

                        // Reddens the models
                        if (this.ebv_MW != 0)
                        {
                            Vector attenuations = unred(Constants.wavelength, -this.ebv_MW); // ebv = 0f - ebv_mw

                            try
                            {
                                model_flux.AddRange((flux * attenuations).Pull());

                            }
                            catch (Exception)
                            {
                                throw new Exception("Spectral Model Run Out of Memory In Fit_Data_To_Models - Get_Models Function");
                            }

                        }
                        else
                        {
                            model_flux.AddRange(flux.Pull());
                        }
                        
                        age_model.Add(Constants.t[i]);
                        metal_model.Add(Constants.Z[j]); // In FireFly 10^{Z[j]} is used to get Metalicity in units of Solar Metalicity

                    }
                }

                this.Model_flux = new Vector(gpu, model_flux.ToArray(), Constants.wavelength.Length); // 2D flattened Array
                this.Model_ages = age_model.ToArray();
                this.Model_metals = metal_model.ToArray();

                return;
            }


        }

        private void TrimModel()
        {

            int length_Data = this.Restframe_Wavelength.Value.Length;
            int length_Mod  = this.Model_wavelength.Length;
            int[] indxs;

            // If data lies within the model
            if (this.Model_wavelength.Value[0] < this.Restframe_Wavelength.Value[0] && this.Model_wavelength.Value[^1] > this.Restframe_Wavelength.Value[^1])
            {

                float[] endVal = (from mwl in this.Model_wavelength.Value
                                  select Math.Abs(mwl - this.Restframe_Wavelength.Value[^1])).ToArray();

                int endIdx = Array.IndexOf(endVal, endVal.Min()) + 1;


                float[] startVal = (from mwl in this.Model_wavelength.Value[0..(endIdx+1)]
                                    select Math.Abs(mwl - this.Restframe_Wavelength.Value[0])).ToArray();

                int startIdx = Array.IndexOf(startVal, startVal.Min());


                if (endIdx - startIdx == length_Data)
                {
                    this.MatchingWavelengthMapping = Enumerable.Range(startIdx, endIdx - startIdx).ToArray();
                    return;
                }


                if (endIdx - startIdx != length_Data)
                {
                    indxs = new int[length_Data];
                    indxs[0] = startIdx;
                    indxs[^1] = endIdx;
                    for (int i = 1; i < indxs.Length-1; i++)
                    {
                        float[] matchtest = (from mwl in this.Model_wavelength.Value[0..endIdx]
                                             select Math.Abs(mwl - this.Restframe_Wavelength[i])).ToArray();
                        indxs[i] = Array.IndexOf(matchtest, matchtest.Min());
                    }

                    this.MatchingWavelengthMapping = indxs;
                    return;
                }

            }

            // If the data lies outside the model
            this.DataLiesOutsideModel = true;

            indxs = new int[length_Data];
            for (int i = 0; i < this.Restframe_Wavelength.Length; i++)
            {
                float[] matchtest = (from mwl in this.Model_wavelength.Value
                                     select Math.Abs(mwl - this.Restframe_Wavelength[i])).ToArray();
                indxs[i] = Array.IndexOf(matchtest, matchtest.Min());
            }
            this.MatchingWavelengthMapping = indxs;

            //Console.WriteLine("WARNING TRIMMODEL FUNC REACHED UNFINISHED CODE");

            return;
        }

        private (int[], Spectral_resolution, VariableGaussianKernel) DownGradeModelInvar(Vector mod_wavelength, int vdisp_round, Vector sres, Vector r_instrument, float min_sig_pix = 0f, bool no_offset = true)
        {
            Vector new_sig = new Vector(gpu, new float[mod_wavelength.Length], 1, false);

            for (int i = 0, j = 0; i < mod_wavelength.Length; i++)
            {
                if (this.MatchingWavelengthMapping[..^1].Contains(i))
                {
                    j++;
                }

                float sig_instrument = Constants.c_kms / (r_instrument.Value[j] * Constants.sig2FWHM);
                new_sig.Value[i] = MathF.Sqrt(MathF.Pow(vdisp_round, 2f) + MathF.Pow(sig_instrument, 2f));

            }
            new_sig.UpdateCache();

            Vector new_sres = Constants.c_div_sig2 / new_sig;

            // IF mod_wavelength < 5 -> raise an error ?? Unlikely

            bool log_wave = true;
            if ((mod_wavelength.Value[3] - mod_wavelength.Value[2]) - (mod_wavelength.Value[2] - mod_wavelength.Value[1]) < 0.000001f * (mod_wavelength.Value[2] - mod_wavelength.Value[1]))
            {
                log_wave = false;
            }


            // match_spectral_resolution Begins

            Vector wave = mod_wavelength;
            Vector new_sres_wave = mod_wavelength;

            if (wave.Min() < new_sres_wave.Value[0] || wave.Max() > new_sres_wave.Value[^1])
            {
                Console.WriteLine("WARNING : Mapping to the new spectral resolution will require extrapolating the provided input vectors!");
            }

            Spectral_resolution new_res = new Spectral_resolution(gpu, new_sres_wave, new_sres, log_wave);
            Spectral_resolution res = new Spectral_resolution(gpu, wave, sres, log_wave);

            res.match(new_res, no_offset, min_sig_pix);



            bool[] mask = Enumerable.Repeat<bool>(false, res.sig_mask.Length).ToArray();

            res.sig_pd.SyncCPU();
            int[] indx = (from sigpd in res.sig_pd.Value
                          where sigpd > min_sig_pix
                          select Array.IndexOf(res.sig_pd.Value, sigpd)).ToArray();




            int len = indx.Length;
            float[] adjusted = res.adjusted_resolution(indx);
            VariableGaussianKernel vGK;


            if (len > 0f)
            {

                float[] selectedsig = new float[len];
                for (int i = 0; i < len; i++)
                {
                    sres.Value[indx[i]] = adjusted[i]; // sres_out
                    selectedsig[i] = res.sig_pd[indx[i]];
                }

                vGK = new VariableGaussianKernel(gpu, selectedsig);

            }
            else
            {

                sres.UpdateCache(adjusted);              // sres_out
                vGK = new VariableGaussianKernel(gpu, res.sig_pd.Pull());
            }

            //if (res.sig_mask.Contains(true) || mask.Contains(true))
            //{
            //    for (int i = 0; i < mask.Length; i++)
            //    {
            //        mask[i] = res.sig_mask[i] || mask[i];  // mask_out
            //    }
            //}



            return (indx, res, vGK);
        }

        private Vector DownGrade(GPU gpu, Vector flux, int[] indx, Spectral_resolution res, VariableGaussianKernel vGK)
        {

            match_spectral_resolution(gpu, flux, indx, res, vGK);
            
            // Call match_spectral_resolution(mod_wavelength, flux, sres, mod_wavelength, new_sres, min_sig_pix=0.0,
            // log10=log_wave, new_log10=log_wave)
            // get : new_flux, matched_sres, sigma_offset, new_mask

            return flux;
        }

        private void match_spectral_resolution(GPU gpu, Vector flux, int[] indx, Spectral_resolution res, VariableGaussianKernel vGK)
        {
            //var dims = 1;

            // 0.2s
            //Spectral_resolution new_res = new Spectral_resolution(gpu, new_sres_wave, new_sres.Value, log10 = new_log10);

            // 0.2s
            //Spectral_resolution res = new Spectral_resolution(gpu, wave, sres, log10);
            // 18-20s
            //res.match(new_res, no_offset, min_sig_pix);
            //float sigma_offset = res.sig_vo;

            // FLUX DEPENDENT SECTION
            //float[] out_flux = flux;

            int len = indx.Length;
            if (len > 0)
            {
                float[] selectedflux = new float[len];
                for (int i = 0; i < len; i++)
                {
                    selectedflux[i] = flux.Value[indx[i]];
                }

                float[] selected_outflux = vGK.Convolve(gpu, new Vector(gpu, selectedflux)).Pull();

                for (int i = 0; i < len; i++)
                {
                    flux.Value[indx[i]] = selected_outflux[i];   // flux out
                }
            }
            else
            {
                //Console.WriteLine($"the length of y is : {flux.Length}");
                flux = vGK.Convolve(gpu, flux);               // flux out
            }
            #region
            // thus ignor else                                      #02

            //float[] adjusted = res.adjusted_resolution(indx);
            //if (len > 0f)
            //{
            //    for (int i = 0; i < indx.Length; i++)
            //    {
            //        out_sres[indx[i]] = adjusted[i];
            //    }
            //}
            //else
            //{
            //    out_sres = adjusted;
            //}

            //if (res.sig_mask.Contains(true) || mask.Contains(true))
            //{
            //    for (int i = 0; i < mask.Length; i++)
            //    {
            //        out_mask[i] = res.sig_mask[i] || mask[i];
            //    }
            //}
            //else
            //{
            //    out_mask = mask;
            //}
            //thus ignore else                                       #01
            #endregion

        }



        private void Match_data_models()
        {
            // Get model_flux and Data_flux




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


        private void NormaliseSpec()
        {
            float data_norm = UtilityMethods.Median(this.Flux.Value);                 
            int num_mods = this.Model_flux.Value.Length / this.Model_flux.Columns;     // 198
            
            Vector model_norm = new Vector(gpu, new float[num_mods],1,false);                                 
            float[] mass_factor = new float[num_mods];                                 


            for (int i = 0; i < num_mods; i++)
            {
                model_norm[i] = UtilityMethods.Median(this.Model_flux.Value[(i*this.Model_flux.Columns)..((i+1) * this.Model_flux.Columns)]); 
                mass_factor[i] = data_norm / model_norm[i];                            
            }

            model_norm.UpdateCache();

            // OVER-WRITES MODEL FLUX WITH THE NORMALISED MODEL FLUX
            this.Model_flux = Vector.ConsecutiveCompoundOperation2D(gpu, this.Model_flux, model_norm, data_norm, "**");
            this.Mass_factor = mass_factor;

        }

        private Vector unred(Vector wavelength, float ebv_mw)
        {




            return wavelength;
        }


        #endregion


        // Test CODE
        public float CalculateChiSqu(int model)
        {
            // Get length of Data
            var length = this.Restframe_Wavelength.Value.Length;

            this.Restframe_Wavelength.SyncCPU();

            // Finds the closest start value
            //float Closest_Start_Val = this.Model_wavelength.OrderBy(n => Math.Abs(this.Restframe_Wavelength.Value[0] - n)).First();

            float Closest_Start_Val = (this.Restframe_Wavelength[0] - this.Model_wavelength).AbsX_IP().Min();


            // Gets Index of Closest Value
            int idx_closest = Array.IndexOf(this.Model_wavelength.Pull(), Closest_Start_Val);

            // Select the Model Flux for the Model
            float[] model_flux = this.Model_flux.Value[((model * this.Model_flux.Columns) + idx_closest)..((model * this.Model_flux.Columns) + idx_closest + length)];

            // Get the Length
            //int length = this.Model_flux.Columns;
            //float[] model_flux = this.Model_flux.Value[(model * length)..((model + 1) * length)];

            float sum = 0f;
            for (int j = 0; j < model_flux.Length; j++)
            {
                sum += MathF.Pow(((model_flux[j] - this.Flux.Value[j]) / this.Error.Value[j]), 2f);
            }

            return sum;
        }
        
        // will combine together SSP 50:50, 30:30:30, 25:25:25:25, 20:20:20:20:20
        private void GenerateInterpSSP()
        {
            int totlen = this.Model_flux.Value.Length;
            int lenOneModel = this.Model_flux.Columns;
            int noModels = totlen / lenOneModel;

            //noModels = 6;

            var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<int>, int>(TrimDownFlux);

            Vector trimmedflux = new Vector(gpu, new float[this.MatchingWavelengthMapping.Length * noModels], this.MatchingWavelengthMapping.Length);
            trimmedflux.IncrementLiveCount();


            MemoryBuffer<float> buffer = trimmedflux.GetBuffer();         // Output
            MemoryBuffer<float> buffer2 = this.Model_flux.GetBuffer();    //  Input
            MemoryBuffer<int> buffer3 = gpu.accelerator.Allocate<int>(this.MatchingWavelengthMapping.Length); //  Input

            buffer3.CopyFrom(gpu.accelerator.DefaultStream, this.MatchingWavelengthMapping, 0, 0, this.MatchingWavelengthMapping.Length);

            kernelWithStream(gpu.accelerator.DefaultStream, noModels, buffer.View, buffer2.View, buffer3.View, lenOneModel);

            gpu.accelerator.Synchronize();

            buffer3.Dispose();

            // RUN 1

            var kernelWithStream2 = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int>(ChiCalcRun1);

            Vector Chis = new Vector(gpu, new float[noModels]);
            Chis.IncrementLiveCount();

            MemoryBuffer<float> buffer_chis = Chis.GetBuffer(); //  Output
            MemoryBuffer<float> buffer_data_flux = this.Flux.GetBuffer(); //  Input
            MemoryBuffer<float> buffer_data_err = this.Error.GetBuffer(); //  Input


            kernelWithStream2(gpu.accelerator.DefaultStream, noModels, buffer_chis.View, buffer_data_flux.View, buffer_data_err.View, buffer.View, this.Flux.Value.Length);

            gpu.accelerator.Synchronize();


            // RUN 2

            List<int> pairX = new List<int>();
            List<int> pairY = new List<int>();

            //int NoG1SSPs = Enumerable.Range(1, noModels).Sum();
            for (int j = noModels - 1, i = 0; j >= 1; j--, i++)
            {
                pairX.AddRange(Enumerable.Repeat(i, j));
                pairY.AddRange(Enumerable.Range(i + 1, j));
            }

            var kernelWithStream3 = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>, ArrayView<int>, int>(ChiCalcRun2);

            Vector ChisG1 = new Vector(gpu, new float[pairX.Count]);
            ChisG1.IncrementLiveCount();

            MemoryBuffer<float> buffer_chis2 = ChisG1.GetBuffer(); //  Output
            MemoryBuffer<int> buffer_X = gpu.accelerator.Allocate<int>(pairX.Count); //  Input
            MemoryBuffer<int> buffer_Y = gpu.accelerator.Allocate<int>(pairY.Count); //  Input

            buffer_chis2.MemSetToZero();

            buffer_X.CopyFrom(gpu.accelerator.DefaultStream, pairX.ToArray(), 0, 0, pairX.Count);
            buffer_Y.CopyFrom(gpu.accelerator.DefaultStream, pairY.ToArray(), 0, 0, pairY.Count);

            kernelWithStream3(gpu.accelerator.DefaultStream, pairX.Count, buffer_chis2.View, buffer_data_flux.View, buffer_data_err.View, buffer.View, buffer_X.View, buffer_Y.View, this.Flux.Value.Length);

            gpu.accelerator.Synchronize();

            buffer_X.Dispose();
            buffer_Y.Dispose();

            //Console.WriteLine($"best model is {Array.IndexOf(Chis, Chis.Min())} with chi : {Chis.Min() / this.Flux.Value.Length}");

            int best_model = Array.IndexOf(ChisG1.Pull(), ChisG1.Min());

            //Console.WriteLine($"\n{this.Path[^25..].Replace(".fits","")} => SSP : {Chis.Min()/this.Flux.Value.Length} || 2 SSP : {ChisG1.Min()/this.Flux.Value.Length}");

            //Console.WriteLine($"best model is {pairX[best_model]} : {pairY[best_model]} with a chi of : {ChisG1.Min() / this.Flux.Value.Length}, with {this.Flux.Value.Length} degrees of freedom");

            RecordData(Array.IndexOf(Chis.Pull(), Chis.Min()), Chis.Min()/this.Flux.Value.Length, new int[2] {pairX[best_model], pairY[best_model]}, ChisG1.Min() / this.Flux.Value.Length);
        } 

        // Kernels
        static void TrimDownFlux(Index1 index, ArrayView<float> Output, ArrayView<float> Fluxes, ArrayView<int> idxs, int fluxlen)
        {
            for (int i = 0; i < idxs.Length; i++)
            {
                Output[index * idxs.Length + i] = Fluxes[index * fluxlen + idxs[i]];
            }
        }
        static void ChiCalcRun1(Index1 index, ArrayView<float> Chis, ArrayView<float> DataFlux, ArrayView<float> DataError, ArrayView<float> ModelFlux, int Modellen)
        {
            float sum = 0f;
            for (int i = 0; i < Modellen; i++)
            {
                sum += XMath.Pow( (ModelFlux[index*Modellen + i] - DataFlux[i]) / DataError[i], 2f);
            }

            Chis[index] = sum;

        }
        static void ChiCalcRun2(Index1 index, ArrayView<float> Chis, ArrayView<float> DataFlux, ArrayView<float> DataError, ArrayView<float> ModelFlux, ArrayView<int> X, ArrayView<int> Y, int Modellen)
        {
            float sum = 0f;
            for (int i = 0; i < Modellen; i++)
            {
                sum += XMath.Pow(( ((ModelFlux[Modellen * X[index] + i] + ModelFlux[Modellen * Y[index] + i]) * 0.5f) - DataFlux[i]) / DataError[i], 2f);
            }

            Chis[index] = sum;
        }


        // record the data output
        private void RecordData(int SSPmodelnum, float SSPchi, int[] dualsspmodelnum, float dualsspchi)
        {
            string data = this.Path.Replace("/Data", "/Output").Replace(".fits","") + ".txt";
            File.WriteAllText(data, "");

            File.AppendAllText(data, "Restframe Wavelength of Data :\n");
            File.AppendAllText(data, string.Join(" ", this.Restframe_Wavelength.Value));

            File.AppendAllText(data, "\nFlux of Data : \n");
            File.AppendAllText(data, string.Join(" ", this.Flux.Value));

            File.AppendAllText(data, "\nError of Data : \n");
            File.AppendAllText(data, string.Join(" ", this.Error.Value));

            File.AppendAllText(data, "\nModel Wavelengths : \n");
            File.AppendAllText(data, string.Join(" ", this.Model_wavelength));


            File.AppendAllText(data, "\n\nBest Single Model Fit : \n");
            File.AppendAllText(data, $"\nModel no : {SSPmodelnum} \n");
            File.AppendAllText(data, $"\nModel chi : {SSPchi} \n");
            File.AppendAllText(data, $"\nGalaxy Age Predicted : {this.Model_ages[SSPmodelnum]} Gyr \n");
            File.AppendAllText(data, $"\nGalaxy Metalicity Predicted : {this.Model_metals[SSPmodelnum]} [Z/H] \n");
            File.AppendAllText(data, $"\nGalaxy Mass Predicted : {(1f / this.Mass_factor[SSPmodelnum]) * 1e-17f} Msolar \n");

            File.AppendAllText(data, "\nModel Flux : \n");
            File.AppendAllText(data, string.Join(" ", this.Model_flux.Value[(SSPmodelnum*this.Model_wavelength.Length)..((SSPmodelnum +1) * this.Model_wavelength.Length)]));



            File.AppendAllText(data, "\n\nBest Dual-SSP Model Fit : \n");
            File.AppendAllText(data, $"\nModel no : {dualsspmodelnum[0]} & {dualsspmodelnum[1]} \n");
            File.AppendAllText(data, $"\nModel chi : {dualsspchi} \n");
            File.AppendAllText(data, $"\nGalaxy Age Predicted : 50% - {this.Model_ages[dualsspmodelnum[0]]} / 50% - {this.Model_ages[dualsspmodelnum[1]]} Gyr \n");
            File.AppendAllText(data, $"\nGalaxy Metalicity Predicted : {0.5f*(this.Model_metals[dualsspmodelnum[0]] + this.Model_metals[dualsspmodelnum[1]])} [Z/H] \n");
            File.AppendAllText(data, $"\nGalaxy Mass Predicted : {0.5f * (((1f / this.Mass_factor[dualsspmodelnum[0]]) * 1e-17f) + ((1f / this.Mass_factor[dualsspmodelnum[1]]) * 1e-17f))} Msolar \n");

            float[] midptflux = Midpoint(
                this.Model_flux.Value[(dualsspmodelnum[0] * this.Model_wavelength.Length)..((dualsspmodelnum[0] + 1) * this.Model_wavelength.Length)],
                this.Model_flux.Value[(dualsspmodelnum[1] * this.Model_wavelength.Length)..((dualsspmodelnum[1] + 1) * this.Model_wavelength.Length)]);

            File.AppendAllText(data, "\nModel Flux : \n");
            File.AppendAllText(data, string.Join(" ", midptflux));



        }

        private float[] Midpoint(float[] vectA, float[] vectB)
        {
            float[] vectR = new float[vectA.Length];
            for (int i = 0; i < vectA.Length; i++)
            {
                vectR[i] = 0.5f * (vectA[i] + vectB[i]);
            }
            return vectR;
        }



        // EXPERIMENTAL CODE DONT USE YET
        public float[] CalculateChiSquVec()
        {
            int length_Data = this.Restframe_Wavelength.Value.Length;
            int length_Mod = this.Model_wavelength.Length;
            int models = length_Data / length_Mod;

            float[] endVal = (from mwl in this.Model_wavelength.Value
                              select Math.Abs(mwl - this.Restframe_Wavelength.Value[^1])).ToArray();

            int endIdx = Array.IndexOf(endVal, endVal.Min()) + 1;


            float[] startVal = (from mwl in this.Model_wavelength.Value[0..(endIdx + 1)]
                                select Math.Abs(mwl - this.Restframe_Wavelength.Value[0])).ToArray();

            int startIdx = Array.IndexOf(startVal, startVal.Min());

            int model_len = endIdx - startIdx;

            List<float> new_mod_flux = new List<float>();
            for (int i = 0; i < models; i++)
            {
                float[] vals = this.Model_flux.Value[(i * length_Mod + startIdx)..(i * length_Mod + endIdx)];
                new_mod_flux.AddRange(vals);
            }


            var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int>(CalcChiKernel);

            Vector output = new Vector(gpu, new float[new_mod_flux.Count]);
            output.IncrementLiveCount();

            Vector new_mod_flux_vec = new Vector(gpu, new_mod_flux.ToArray());
            new_mod_flux_vec.IncrementLiveCount();

            MemoryBuffer<float> buffer = output.GetBuffer(); // Output
            MemoryBuffer<float> buffer2 = new_mod_flux_vec.GetBuffer(); // Input
            MemoryBuffer<float> fluxBuf = Flux.GetBuffer(); // Input
            MemoryBuffer<float> errorBuf = Error.GetBuffer(); // Input


            kernelWithStream(gpu.accelerator.DefaultStream, buffer.Length, buffer.View, buffer2.View, fluxBuf.View, errorBuf.View , model_len);

            gpu.accelerator.Synchronize();

            output.SyncCPU();

            float[] chis = new float[models];
            for (int i = 0; i < models; i++)
            {
                chis[i] = output.GetRowAsArray(i,true).Sum();
            }

            return chis;
        }
        // Kernel
        static void CalcChiKernel(Index1 index, ArrayView<float> Output, ArrayView<float> InputMod_flux, ArrayView<float> InputDat_flux, ArrayView<float> InputDat_err, int model_len)
        {
            Output[index] = XMath.Pow((InputMod_flux[index] - InputDat_flux[index % model_len]) / InputDat_err[index % model_len], 2f);
        }




    }



}
