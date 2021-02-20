using ILGPU;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;



namespace MachineLearningSpectralFittingCode
{
    class Program : UtilityMethods
    {
        // Singleton Instances
        public static Config config;
        public static Cosmology cosmology;


        static void Main()
        {

            // PRE-INITIALISATION
            config = new Config();
            config.GetHardware();
            cosmology = new Cosmology();
            cosmology.Initialise();


            //var watch = System.Diagnostics.Stopwatch.StartNew();



            // VARIABLE BLOCK
            string Data_path = @"C:\Users\Marcelpaw\source\repos\MachineLearningSpectralFittingCode\Data\spec-0266-51602-0001.dat";

            // PROGRAM START
            Console.WriteLine("Start");






            // READ IN DATA
            //Vector Data = new Vector(ReadData(Data_path), 3); // Data Is read in as a 2D Vector of 3 columns

            // Spectra Spectra = new Spectra(Data_path, Config.Milky_Way_Reddening, Config.HPF_Mode, Config.N_Masked_Amstrongs);
            //Spectra = Spectra.InitialiseSpectraParameters(AcceleratorIds[1], Spectra, Data, Config.Redshift, Config.RA_DEC, Config.Velocity_Dispersion, Config.Instrument_Resolution);


            //var watch = System.Diagnostics.Stopwatch.StartNew();
            float ans = cosmology.luminosity_distance(0.99f);
            Console.WriteLine($"the luminosity distance at 0.99 redshift is {ans} ");
            //watch.Stop();
            //var elapsedMs = watch.ElapsedMilliseconds;
            //Console.WriteLine("Time Taken to reach setup " + (elapsedMs * 0.001f).ToString() + "s");
            //Console.WriteLine($"the luminosity distance at 0.99 redshift is {ans} ");

            //for (int i = 0; i < TestVectorB.Value.Length; i++)
            //{
            //    Console.WriteLine(TestVectorB.Value[i]);
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

        }



    }
}
