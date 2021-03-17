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
        public static string PathOfProgram = "C:/Users/Marcelpaw/source/repos/MachineLearningSpectralFittingCode/";
        public static Context context;
        public static Accelerator gpu;

        static void Main()
        {
            var watch = System.Diagnostics.Stopwatch.StartNew();

            // PRE-INITIALISATION
            File.WriteAllText($"{PathOfProgram}Log.txt", $"{System.DateTime.Now} : Starting Pre-Initialisation\n");
            config = new Config();
            cosmology = new Cosmology();
            cosmology.Initialise();
            //config.RecordSystemInfo();





            // VARIABLE BLOCK
            string Data_path = PathOfProgram + @"\Data\spec-0266-51602-0001.dat";

            // PROGRAM START
            Console.WriteLine("Start");


            
            //var tree = Hdf5.ReadTreeFileStructure(fileName);

            //float[] fluxgrid = new float[4563];
            //for (int i = 0; i < 4563; i++)
            //{
            //    fluxgrid[i] = readobj[0, 0, 0, i];
            //}

            //Console.WriteLine(fluxgrid.ToString());


            Console.WriteLine("reading hdf5");
            //for (int i = 0; i < hdata.Length; i++)
            //{
            //    Console.WriteLine(hdata[i]);
            //}



            //// READ IN DATA
            Vector Data = new Vector(UtilityMethods.ReadData(Data_path), 3); // Data Is read in as a 2D Vector of 3 columns
            



            Spectral_Model spectral_Model = new Spectral_Model(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs);
            spectral_Model.InitialiseSpectraParameters(gpu, Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
            spectral_Model.Fit_models_to_data();

            float[] chis = new float[spectral_Model.Model_ages.Length];

            for (int i = 0; i < chis.Length; i++)
            {
                chis[i] = spectral_Model.CalculateChiSqu(i);
            }



            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;
            Console.WriteLine("Time Taken to reach setup " + (elapsedMs * 0.001f).ToString() + "s");

            for (int i = 0; i < chis.Length; i++)
            {
                Console.WriteLine($"model {i} : chi {chis[i]}");
            }


            // Get length of Data
            var length = spectral_Model.Restframe_Wavelength.Value.Length;

            // Finds the closest start value
            float Closest_Start_Val = spectral_Model.Model_wavelength.OrderBy(n => Math.Abs(spectral_Model.Restframe_Wavelength.Value[0] - n)).First();

            // Gets Index of Closest Value
            int idx_closest = Array.IndexOf(spectral_Model.Model_wavelength, Closest_Start_Val);



            // WRITES MODEL DATA TO HDF5 FORMAT FOR PLOTTING
            string fileName = Program.PathOfProgram + @"modelOutput.h5";
            long fileId = Hdf5.CreateFile(fileName);

            Hdf5.WriteDatasetFromArray<float>(fileId, "model_wavelength_trimmed", spectral_Model.Model_wavelength[idx_closest..(idx_closest+length)]);
            
            Hdf5.WriteDatasetFromArray<float>(fileId, "model_flux_norm", spectral_Model.Model_flux.Value);
            Hdf5.WriteDatasetFromArray<int>(fileId, "model_flux_shape", new int[3] { spectral_Model.Model_ages.Length, spectral_Model.Model_flux.Columns, idx_closest });
            
            Hdf5.WriteDatasetFromArray<float>(fileId, "model_ages", spectral_Model.Model_ages);
            Hdf5.WriteDatasetFromArray<float>(fileId, "model_metals", spectral_Model.Model_metals);
            Hdf5.WriteDatasetFromArray<float>(fileId, "model_chis", chis);

            Hdf5.WriteDatasetFromArray<float>(fileId, "data_wavelength_RestFrame", spectral_Model.Restframe_Wavelength.Value);
            Hdf5.WriteDatasetFromArray<float>(fileId, "data_flux", spectral_Model.Flux.Value);
            Hdf5.WriteDatasetFromArray<float>(fileId, "data_error", spectral_Model.Error.Value);

            Hdf5.CloseFile(fileId);


            // Read Models is 


            // MW & Background Reddening
            // Dust attenuation


            // UI

            // LEARNING

            // ANALYSIS

            // OUTPUT

            // RUN PYTHON FOR VISUALS

            // PROGRAM END
            Console.WriteLine("End");
            File.AppendAllText($"{PathOfProgram}Log.txt", $"{System.DateTime.Now} : Program Terminated \n");

        }



    }
}
