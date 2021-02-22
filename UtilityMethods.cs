using System;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;


namespace MachineLearningSpectralFittingCode
{
    public class UtilityMethods
    {

        public static float[] ReadData(string Path)
        {

            string readText = File.ReadAllText(Path);

            readText = readText.Replace("\n", " ");
            float[] Data = Array.ConvertAll(readText.Split(' '), float.Parse);

            //for (int i = 0; i < Data.Length/2; i++)
            //{
            //    Console.WriteLine(Data[i]);
            //}

            return Data;
        } 

        public static float Degree2Radians(float Deg)
        {
            return Deg * Constants.Deg2RadFactor;
        }

        public static float Mpc2cm(float Mpc)
        {
            return 3.08567758128e+24f * Mpc;
        }


    }
}
