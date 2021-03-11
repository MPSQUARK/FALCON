using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    public class SPModel
    {
        int
            velocity_dispersion_r
        { get; set; }
        int
            fit_per_iteration_cap
        { get; set; }

        List<double>
            delta_lamdba_lib
        { get; set; }

        public void InitialiseSPModel(float velocity_disp)
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
            this.velocity_dispersion_r = (int)(MathF.Round(velocity_disp / 5f) * 5f);

        }

        public void fit_models_to_data_Calc(float velocity_dispersion, Vector restframe_wavelength, Vector r_instrument, float ebv_MW)
        {
            // enumerate over each model_lib (mi, mm)?
            // enumerate over imfs (ii)?
            //var deltal = this.delta_lamdba_lib[mi];
            Get_Model(Program.config.Model_Key, Program.config.IMF, this.delta_lamdba_lib[0], velocity_dispersion, restframe_wavelength, r_instrument, ebv_MW);
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

        private void Get_Model(byte Model_Key, byte IMF, double delta_l, float velocity_dispersion, Vector restframe_wavelength, Vector r_instrument, float ebv_MW)
        {
            if (Model_Key % 2 == 1)
            {
                var first_file = true;
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
            else if (Model_Key % 2 == 0)
            {

                var model_path = Program.PathOfProgram + @"\StellarPopulationModels\MaStar_SSP_v0.2.fits";

                // lib = model_used
                var slope = 0f;
                if (IMF == 0) // kr
                {
                    slope = 1.3f;
                }
                else if (IMF == 1) // ss
                {
                    slope = 2.35f;
                }
                else
                {
                    throw new Exception("Unrecognised IMF");
                }

                var sidx = Array.IndexOf(Constants.s, slope);
                
                Console.WriteLine(sidx);
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

    }


}
