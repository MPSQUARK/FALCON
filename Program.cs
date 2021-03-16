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
            //var watch = System.Diagnostics.Stopwatch.StartNew();

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
            

            //// Made 1 Instance of a Spectrum

            //Parallel.For(0, 1, i =>
            //{
            //Spectra Spectrum = new Spectra(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs);
            //Spectrum.InitialiseSpectraParameters(gpu, Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
            //Spectrum.Fit_models_to_data();
            //    Spectrum = null;
            //}
            //);

            //Spectral_Model spectral_Model = new Spectral_Model(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs);
            //spectral_Model.InitialiseSpectraParameters(gpu, Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
            //spectral_Model.Fit_models_to_data();

            Console.WriteLine();

            Vector testvec = new Vector(new float[20] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 }, 5);
            Vector testvec2 = new Vector(new float[5] { 5, 10, 15, 20, 25 }, 1);
            Vector ans = Vector.ConsecutiveOperation2D(gpu, testvec, testvec2, '*');

            //for (int j = 0; j < ans.Value.Length / ans.Columns; j++)
            //{
            //    for (int i = 0; i < ans.Columns; i++)
            //    {
            //        Console.Write($" {ans.Value[i + j * ans.Columns]} |");
            //    }
            //    Console.WriteLine();
            //}

            for (int i = 0; i < ans.Value.Length; i++)
            {
                Console.WriteLine(ans.Value[i]);
            }


            //var chi = spectral_Model.CalculateChiSqu(0);
            //Console.WriteLine($"The Chi Squared Of Model : {0} , is {chi}");

            //float chi = 0f;
            //chi = spectral_Model.CalculateChiSqu(0);
            //Console.WriteLine(chi);
            //chi = spectral_Model.CalculateChiSqu(150);
            //Console.WriteLine(chi);



            //chi = spectral_Model.CalculateChiSqu(10);
            //chi = spectral_Model.CalculateChiSqu(20);
            //chi = spectral_Model.CalculateChiSqu(50);
            //chi = spectral_Model.CalculateChiSqu(100);

            // WRITES MODEL DATA TO HDF5 FORMAT FOR PLOTTING
            //string fileName = Program.PathOfProgram + @"modeldata.h5";
            //long fileId = Hdf5.CreateFile(fileName);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_wavelength", spectral_Model.Model_wavelength);
            //Hdf5.WriteDatasetFromArray<int>(fileId, "model_flux_shape", new int[2] { spectral_Model.Model_flux.Length, spectral_Model.Model_flux[0].Length });
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_flux", spectral_Model.Model_flux.SelectMany(a => a).ToArray());
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_ages", spectral_Model.Model_ages);
            //Hdf5.WriteDatasetFromArray<float>(fileId, "model_metals", spectral_Model.Model_metals);


            //var slice = Vector.AccessSlice(gpu, Data, 0, 'c');
            //for (int i = 0; i < slice.Value.Length; i++)
            //{
            //    Console.WriteLine(slice.Value[i]);
            //}


            //watch.Stop();
            //var elapsedMs = watch.ElapsedMilliseconds;
            //Console.WriteLine("Time Taken to reach setup " + (elapsedMs * 0.001f).ToString() + "s");


            //Console.WriteLine(Spectrum.Redshift);
            //Console.WriteLine(Spectrum.Distance_Luminosity);

            // Made 2 Instance of a Spectrum
            //Spectra Spectrum2 = new Spectra(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs);
            //Spectrum2.InitialiseSpectraParameters(config.AcceleratorIds[config.GPU_ids[0]], Data, 0.51f, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);

            //Console.WriteLine(Spectrum2.Redshift);
            //Console.WriteLine(Spectrum2.Distance_Luminosity);

            //float ans = cosmology.luminosity_distance(0.99f);
            //Console.WriteLine($"the luminosity distance at {Spectrum.Redshift} redshift is {Spectrum.Distance_Luminosity} ");
            //Console.WriteLine($"the luminosity distance at {0.99f} redshift is {ans} ");

            //Console.WriteLine($"the luminosity distance at 0.99 redshift is {ans} ");

            //for (int i = 0; i < Spectrum.Wavelength.Value.Length; i++)
            //{
            //    Console.WriteLine(Spectrum.Wavelength.Value[i]);
            //}


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
