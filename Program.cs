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
            //string Data_path = @"C:\Users\Marcelpaw\source\repos\MachineLearningSpectralFittingCode\spec-0266-51602-0001.dat";

            //ReadData(Data_path);
            Vector vectorA = new Vector(new float[3] { 1, 3, -5 });
            Vector vectorB = new Vector(new float[3] { 4, -2, -1 });
            float scalar = 5f;

            List<AcceleratorId> AcceleratorIds = new List<AcceleratorId>();

            foreach (var accelerator in Accelerator.Accelerators)
            {
                AcceleratorIds.Add(accelerator);
            }

            //Vector Output = Vector.ScalarProduct(AcceleratorIds[1] ,vectorA, scalar);


            //for (int i = 0; i < Output.value.Length; i++)
            //{
            //    Console.WriteLine(Output.value[i]);
            //}

            float test_val =  Vector.DotProduct(AcceleratorIds[1], vectorA, vectorB);
            Console.WriteLine(test_val);
            
            

            Console.WriteLine("End");






            //Console.WriteLine("Start");


            // Read Data in 
            // Read Modles is 


            // MW & Background Reddening
            // Dust attenuation


            // UI

            // LEARNING

            // ANALYSIS

            // OUTPUT

            // RUN PYTHON FOR VISUALS


        }
    }
}
