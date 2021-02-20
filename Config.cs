using System;
using ILGPU.Runtime;
using System.Collections.Generic;
using System.Text;

using System.Management; // Uninstall later when system analytics unnessessary
using System.IO;// Uninstall later when system analytics unnessessary
using ILGPU;



namespace MachineLearningSpectralFittingCode
{
    public class Config
    {
        // Constructor
        public Config()
        {
            GetHardware();
        }

        // CONFIG OF HARDWARE   
        public List<AcceleratorId> AcceleratorIds = new List<AcceleratorId>();
        public List<byte> GPU_ids = new List<byte>();
        public bool HasGPU = false;
        public bool Has_Multi_GPU  = false;

        public void GetHardware()
        {

            foreach (var accelerator in Accelerator.Accelerators)
            {

                string type = accelerator.AcceleratorType.ToString();
                float id = 0;
                if (type == "Cuda" || type == "OpenCL")
                {
                    this.AcceleratorIds.Add(accelerator);
                    this.GPU_ids.Add((byte)id);
                    id++;
                }
                else if (type == "CPU")
                {
                    continue;
                }
                else
                {
                    Console.WriteLine("Unknown hardware detected");
                }
            }

            if (this.GPU_ids.Count > 1)
            {
                this.Has_Multi_GPU = true;
            }
            else if (GPU_ids.Count == 1)
            {
                this.HasGPU = true;
            }
            else
            {
                Console.WriteLine("Warning no GPU detected");
            }
        }

        // To be used for code performace analysis accross various hardware
        public void RecordSystemInfo()
        {

            // Console.WriteLine("64 Bit operating system? : {0}", Environment.Is64BitOperatingSystem ? "Yes" : "No");

            ManagementClass myManagementClass = new ManagementClass("Win32_Processor");
            ManagementObjectCollection myManagementCollection = myManagementClass.GetInstances();
            PropertyDataCollection myProperties = myManagementClass.Properties;

            foreach (var obj in myManagementCollection)
            {

                File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"Device Name : {obj.Properties["SystemName"].Value}".Trim()+"\n");
                File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"CPU Model : {obj.Properties["Name"].Value}".Trim() +"\n");
                File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"CPU Base Frequency : {obj.Properties["CurrentClockSpeed"].Value}".Trim() +"\n");
                File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"CPU Thread Count : {obj.Properties["ThreadCount"].Value}".Trim()+"\n");

            }

            // Create main context
            using (var context = new Context())
            {
                // Create default accelerator for the given accelerator id.
                // Note that all accelerators have to be disposed before the global context is disposed
                using (var accelerator = Accelerator.Create(context, AcceleratorIds[GPU_ids[0]]))
                {

                    File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"GPU Name : {accelerator.Name}".Trim() + "\n");
                    File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"GPU Memory : {accelerator.MemorySize / MathF.Pow(1024f, 3f)}".Trim() + "\n");
                    File.AppendAllText($"{Program.PathOfProgram}Log.txt", $"GPU Cores : {accelerator.MaxNumThreadsPerGroup}".Trim() + "\n");

                    accelerator.Dispose();
                }
                context.Dispose();
            }


        }



        // CONFIG OF PHYSICS
        #region

        // Cosmology Model Set
        public static Cosmology cosmology { get; private set; } = new Cosmology();


        // Redshift Value
        public float Redshift { get; private set; } = 0.021275453f;
        // Right Ascension and Declination
        public float[] RA_DEC { get; private set; } = new float[2] { 145.89219f, 0.059372f };
        // Velocity Dispersion in km/s
        public float Velocity_Dispersion { get; private set; } = 135.89957f;
        // Instrument Resolution !!errorNES!!
        public float Instrument_Resolution { get; private set; } = 2000f;
        // Number masked amstrongs !!errorNES!!
        public ushort N_Masked_Amstrongs { get; private set; } = 20;
        // Choose Model - DEFAULT = 0 : MaStar, 1 : m11
        public byte Model_Key { get; private set; } = 0;
        // Choose Model Flavour 
        // FOR Model 0 - DEFAULT = 0 : E-MaStar, 1 : Th-MaStar
        // FOR Model 1 - DEFAULT = 0 : MILES, 1 : STELIB, 2 : ELODIE, 3 : MARCS (kr imf only)
        public byte Model_Flavour { get; private set; } = 0;
        // Choose IMF : DEFAULT 0 : Kroupa, 1 : Salpeter
        public byte IMF { get; private set; } = 0;
        // Min and Max Age of Models
        public float[] MinMax_Model_Age { get; private set; } = new float[2] { 0, Constants.AoU };
        // Min and Max Metalicity of Models
        public float[] MinMax_Model_Metalicity { get; private set; } = new float[2] { -3f, 3f };
        // Is the Data in Air or Vacuum - DEFAULT true : Air, false : Vacuum
        public bool Data_Medium { get; private set; } = false;
        // Flux Scale Factor, this code uses units of erg/s/A/cm^2
        // e.g. Flux_Scale_Factor = -17 for SDSS Data which means the data is scaled by 10^(-17)
        // Value range -128 to 127
        public sbyte Flux_Scale_Factor { get; private set; } = -17;
        // Write Results to Output file
        public bool Write_Results { get; private set; } = true;
        // Correct for Milky Way Reddening
        public bool Milky_Way_Reddening { get; private set; } = true;
        // Set Parameters for dust determination - DEFAULT true : 'on', false : 'hpf only' i.e E(B-V)=0
        public bool HPF_Mode { get; private set; } = true;
        // Set Dust Law - DEFAULT 0 : Calzetti, 1 : Allen, 2 : prevot
        public byte Dust_Law { get; private set; } = 0;

        #endregion

        // CONFIG OF AI

    }
}
