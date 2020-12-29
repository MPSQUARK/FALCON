using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using System.Globalization;

namespace MachineLearningSpectralFittingCode
{
    public class ProcessDataMethods
    {

        public static void ReadData(string Path)
        {



            string readText = File.ReadAllText(Path);
            string[] text = readText.Split(' ', '\n');

            float[] result = Array.ConvertAll(text,float.Parse);

            Console.WriteLine(result);

            return;
        } 



    }
}
