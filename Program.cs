using ILGPU;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;


namespace MachineLearningSpectralFittingCode
{
    class Program : ProcessDataMethods
    {

        static void Main()
        {
            // PRE-INITIALISATION VARIABLE BLOCK
            List<AcceleratorId> AcceleratorIds = new List<AcceleratorId>();
            byte IndexOfVectorDevice = 1;

            // PRE-INITIALISATION
            foreach (var accelerator in Accelerator.Accelerators)
            {
                AcceleratorIds.Add(accelerator);
            }

            // VARIABLE BLOCK
            string Data_path = @"C:\Users\Marcelpaw\source\repos\MachineLearningSpectralFittingCode\spec-0266-51602-0001.dat";

            // PROGRAM START
            Console.WriteLine("Start");

            // READ IN DATA
            Vector Data = new Vector(ReadData(Data_path), 3); // Data Is read in as a 2D Vector of 3 columns

            Vector DataCol0 = Vector.AccessSlice(AcceleratorIds[IndexOfVectorDevice], Data, 0, 'c');

            //for (int i = 0; i < DataCol0.Value.Length; i++)
            //{
            //    Console.WriteLine(DataCol0.Value[i]);
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
