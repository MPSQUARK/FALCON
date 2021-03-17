﻿using System; // System Stuff
using System.Collections.Generic;
using System.Text;

using ILGPU;  // GPU MODULE 
using ILGPU.Runtime; // GPU MODULE

using System.Management; // Uninstall later when system analytics unnessessary
using System.IO;// Uninstall later when system analytics unnessessary

using HDF5CSharp;

namespace MachineLearningSpectralFittingCode
{
    public class Config
    {
        // Constructor
        public Config()
        {
            // Get Hardware Data
            GetHardware();
            // Get Model Data 
            GetModelData();
        }

        // CONFIG OF HARDWARE   
        private List<AcceleratorId> AcceleratorIds = new List<AcceleratorId>();
        private List<byte> GPU_ids = new List<byte>();
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

            Program.context = new Context();
            Program.context.EnableAlgorithms();
            Program.gpu = Accelerator.Create(Program.context, AcceleratorIds[GPU_ids[0]]);
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

        /* Model Key Selector 
         * first 3 bits denote the Model type 
         * last  5 bits denote the Model flavour
         * e.g. 0b00000_001 => m11
         * e.g. 0b00001_001 => m11 - Miles
         * e.g. 0b00010_001 => m11 - Stelib
         * e.g. 0b00100_001 => m11 - Elodie
         * e.g. 0b01000_001 => m11 - Marcs (kr imf only)
         * e.g. 0b10000_001 => m11 - **NO ENTRY**
         * e.g. 0b00000_010 => MaStar
         * e.g. 0b00001_010 => MaStar - Th
         * e.g. 0b00010_010 => MaStar - E
         * e.g. 0b00000_100 => **NO ENTRY** 
        */
        public byte Model_Key { get; private set; } = 0b00010_010;

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
        public bool Downgrade_models { get; private set; } = true;

        #endregion

        // CONFIG OF AI


        // CONFIG PRE-INITIALISE DATA
        public void GetModelData()
        {
            // run this function upon Config Application - Pre-Initialisation Phase

            // Read in Hdf5 Data File/s
            string fileName = Program.PathOfProgram + @"data.h5";
            long id = Hdf5.OpenFile(fileName, true);

            // Load Data For Ma-Star Models
            if (Model_Key % 2 == 0)
            {
                try
                {
                    // READ IN DATA
                    Constants.r_model = (double[])((Hdf5.ReadDatasetToArray<double>(id, "r_model")).result);
                    
                    Constants.t = (float[])((Hdf5.ReadDatasetToArray<float>(id, "t")).result);
                    Constants.Z = (float[])((Hdf5.ReadDatasetToArray<float>(id, "Z")).result);
                    Constants.s = (float[])((Hdf5.ReadDatasetToArray<float>(id, "s")).result);

                    Constants.wavelength = (float[])((Hdf5.ReadDatasetToArray<float>(id, "wavelength")).result);

                    // READ IN FLUXGRID
                    float[,,,] fluxgrid;
                    if (Model_Key == 0b00001_010) // MaStar-Th
                    {
                        fluxgrid = (float[,,,])((Hdf5.ReadDatasetToArray<float>(id, "fluxgrid_Th")).result);
                    }
                    else if (Model_Key == 0b00010_010) // MaStar-E
                    {
                        fluxgrid = (float[,,,])((Hdf5.ReadDatasetToArray<float>(id, "fluxgrid_E")).result);
                    }
                    else
                    {
                        throw new Exception("MaStar Model Flavour Error Please choose between MaStar-Th and MaStar-E");
                    }

                    Hdf5.CloseFile(id);


                    // SET FLUXGRID 
                    Constants.fluxgrid = fluxgrid;


                    // SET SLOPE AND SIDX
                    switch (this.IMF)
                    {
                        case 0: // kr
                            Constants.slope = 1.3f;
                            break;
                        case 1: // ss
                            Constants.slope = 2.35f;
                            break;

                        default:
                            throw new Exception("Unrecognised IMF");
                    }
                    Constants.sidx = Array.IndexOf(Constants.s, Constants.slope);

                }
                catch (Exception)
                {
                    throw new Exception("MaStar SSP Data not found");
                }
            }
            
        }

    }
}
