using ILGPU;
using ILGPU.Runtime;
using System;
using System.IO;
using System.Collections.Generic;
using System.Threading.Tasks;
using HDF5CSharp;

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


            //// Made 1 Instance of a Spectrum

            //Parallel.For(0, 1, i =>
            //{
                Spectra Spectrum = new Spectra(Data_path, config.Milky_Way_Reddening, config.HPF_Mode, config.N_Masked_Amstrongs);
                Spectrum.InitialiseSpectraParameters(gpu, Data, config.Redshift, config.RA_DEC, config.Velocity_Dispersion, config.Instrument_Resolution);
                Spectrum.Fit_models_to_data();
            //    Spectrum = null;
            //}
            //);

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
