using System;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;


namespace MachineLearningSpectralFittingCode
{
    public class ProcessDataMethods
    {

        public static float[][] ReadData(string Path)
        {

            string readText = File.ReadAllText(Path);
            // Row-Column format
            float[][] Data = readText.Split('\n')
                .Select(x => x.Split(' '))
                .Select(x => Array.ConvertAll(x, float.Parse))
                .ToArray();


            //for (int j = 0; j < Data.Length; j++)
            //{
            //    for (int i = 0; i < Data[0].Length; i++)
            //    {
            //        Console.Write(Data[j][i]);
            //        Console.Write(" , ");
            //    }
            //    Console.Write("\n");
            //}

            return Data;
        } 



    }
}
