using ILGPU;
using ILGPU.Runtime;
using System;
using System.IO;
using System.Collections.Generic;
using System.Threading.Tasks;
using HDF5CSharp;
using System.Linq;


namespace MachineLearningSpectralFittingCode
{
    class Program
    {
        // Singleton Instances
        public static Config config;
        public static Cosmology cosmology;
        public static string PathOfProgram = AppDomain.CurrentDomain.BaseDirectory; //"C:/Users/Marcelpaw/source/repos/MachineLearningSpectralFittingCode/";

        static void Main()
        {
            // Variable BLOCK
            Context context = new Context();
            context.EnableAlgorithms();
            Accelerator gpu;

            config = new Config();
            gpu = config.GetHardware(context);
            
            cosmology = new Cosmology();
            cosmology.Initialise();

            //UI UserInterface = new UI(config);
            config.GetModelData(gpu);

            string Data_path = PathOfProgram + @"\Data\spec-0266-51602-0001.dat";



            // Timer
            var watch = System.Diagnostics.Stopwatch.StartNew();

            /* PRE-INITIALISATION
            */
            File.WriteAllText($"{PathOfProgram}Log.txt", $"{System.DateTime.Now} : Starting Pre-Initialisation\n");

            // Set Configs


            //config.RecordSystemInfo();

            // Set Cosmology






            // Enable UI to SET User Defined config


            // PROGRAM START
            Console.WriteLine("Start");


            //// READ IN DATA
            Vector Data = new Vector(UtilityMethods.ReadData(Data_path), 3); // Data Is read in as a 2D Vector of 3 columns

            Spectral_Model spectral_Model = new Spectral_Model(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs, gpu);
            spectral_Model.InitialiseSpectraParameters(Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
            spectral_Model.Fit_models_to_data();

            float[] chis = new float[spectral_Model.Model_ages.Length];

            for (int j = 0; j < chis.Length; j++)
            {
                chis[j] = spectral_Model.CalculateChiSqu(j);
                Console.WriteLine(chis[j] / spectral_Model.Flux.Value.Length);
            }

            Console.WriteLine($"best model is {Array.IndexOf(chis, chis.Min())}");




            //Spectral_Model[] spectral_Models = new Spectral_Model[500];

            //Parallel.For(0, 500, i =>
            //{

            //   spectral_Models[i] = new Spectral_Model(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs, gpu);
            //   spectral_Models[i].InitialiseSpectraParameters(Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
            //   spectral_Models[i].Fit_models_to_data();

            //   float[] chis = new float[spectral_Models[i].Model_ages.Length];

            //for (int j = 0; j < chis.Length; j++)
            //{
            //    chis[j] = spectral_Models[i].CalculateChiSqu(j);
            //}

            ////float[] chisFast = spectral_Model.CalculateChiSquVec(gpu);
            //});


            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;
            Console.WriteLine("Time Taken to Complete " + (elapsedMs * 0.001f).ToString() + "s");
            Console.WriteLine("Time Taken to Complete per Spectra " + (elapsedMs * 0.001f / 500f).ToString() + "s");


            //for (int i = 0; i < chis.Length; i++)
            //{
            //    Console.WriteLine($"{chis[i] / spectral_Model.Flux.Value.Length} index {i}"); //- chisFast[i]);
            //}



            //Console.WriteLine("Recording System Info");


            //config.RecordSystemInfo();


            /* GETS THE MASS OF THE GALAXY
             * Console.WriteLine(((1f/spectral_Model.Mass_factor[188]) * 1e-17f).ToString());
             */


            // Get length of Data
            //var length = spectral_Model.Restframe_Wavelength.Value.Length;

            //// Finds the closest start value
            //float Closest_Start_Val = spectral_Model.Model_wavelength.OrderBy(n => Math.Abs(spectral_Model.Restframe_Wavelength.Value[0] - n)).First();

            //// Gets Index of Closest Value
            //int idx_closest = Array.IndexOf(spectral_Model.Model_wavelength, Closest_Start_Val);



            // WRITES MODEL DATA TO HDF5 FORMAT FOR PLOTTING
            //string fileName = Program.PathOfProgram + @"modelOutput.h5";
            //long fileId = Hdf5.CreateFile(fileName);

            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_wavelength_trimmed", spectral_Model.Model_wavelength[idx_closest..(idx_closest+length)]);

            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_flux_norm", spectral_Model.Model_flux.Value);
            //Hdf5.WriteDatasetFromArray<int>(fileId, "model_flux_shape", new int[3] { spectral_Model.Model_ages.Length, spectral_Model.Model_flux.Columns, idx_closest });

            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_ages", spectral_Model.Model_ages);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_metals", spectral_Model.Model_metals);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_chis", chis);

            //Hdf5.WriteDatasetFromArray<float>(fileId, "data_wavelength_RestFrame", spectral_Model.Restframe_Wavelength.Value);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "data_flux", spectral_Model.Flux.Value);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "data_error", spectral_Model.Error.Value);

            //Hdf5.CloseFile(fileId);


            // Read Models is 


            // MW & Background Reddening
            // Dust attenuation


            // UI

            // LEARNING

            // ANALYSIS

            // OUTPUT

            // RUN PYTHON FOR VISUALS

            // PROGRAM END
            //Console.WriteLine("End");
            //File.AppendAllText($"{PathOfProgram}Log.txt", $"{System.DateTime.Now} : Program Terminated \n");

            //Console.ReadLine();

        }



    }
}
