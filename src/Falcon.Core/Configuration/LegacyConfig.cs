using HDF5CSharp;
using BAVCL;
using Falcon.Core.Constants;
using Falcon.Core.Services.Cosmology;

namespace Falcon.Core.Configuration;

public static class LegacyConfig
{
    public static string PathOfProgram = AppDomain.CurrentDomain.BaseDirectory;
    public static CosmologyService Cosmology = new();

    // CONFIG OF HARDWARE   
    public static bool Has_Multi_GPU = false;

    // CONFIG OF PHYSICS
    #region


    // Redshift Value
    public static float Redshift { get; private set; } = 0.021275453f;
    // Right Ascension and Declination
    public static float[] RA_DEC { get; private set; } = [145.89219f, 0.059372f];
    // Velocity Dispersion in km/s
    public static float Velocity_Dispersion { get; private set; } = 135.89957f;
    // Instrument Resolution !!errorNES!!
    public static float Instrument_Resolution { get; private set; } = 2000f;
    // Number masked amstrongs !!errorNES!!
    public static ushort N_Masked_Amstrongs { get; private set; } = 20;

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
    public static byte Model_Key = 0b00010_010;

    // Choose IMF : DEFAULT 0 : Kroupa, 1 : Salpeter
    public static byte IMF { get; private set; } = 0;
    // Min and Max Age of Models
    public static float[] MinMax_Model_Age { get; private set; } = [0, FireflyConstants.AoU];
    // Min and Max Metalicity of Models
    public static float[] MinMax_Model_Metalicity { get; private set; } = [-3f, 3f];
    // Is the Data in Air or Vacuum - DEFAULT true : Air, false : Vacuum
    public static bool Data_Medium { get; private set; } = false;
    // Flux Scale Factor, this code uses units of erg/s/A/cm^2
    // e.g. Flux_Scale_Factor = -17 for SDSS Data which means the data is scaled by 10^(-17)
    // Value range -128 to 127
    public static sbyte Flux_Scale_Factor { get; private set; } = -17;
    // Write Results to Output file
    public static bool Write_Results { get; private set; } = true;
    // Correct for Milky Way Reddening
    public static bool Milky_Way_Reddening { get; private set; } = true;
    // Set Parameters for dust determination - DEFAULT true : 'on', false : 'hpf only' i.e E(B-V)=0
    public static bool HPF_Mode { get; private set; } = true;
    // Set Dust Law - DEFAULT 0 : Calzetti, 1 : Allen, 2 : prevot
    public static byte Dust_Law { get; private set; } = 0;
    public static bool Downgrade_models { get; private set; } = true;

    #endregion

    public static void ApplyFrom(FireflyOptions options, string basePath)
    {
        PathOfProgram = basePath.EndsWith(Path.DirectorySeparatorChar) ? basePath : basePath + Path.DirectorySeparatorChar;
        Model_Key = options.Model.ModelKey;
        IMF = options.Model.Imf.Equals("Salpeter", StringComparison.OrdinalIgnoreCase) ? (byte)1 : (byte)0;
        MinMax_Model_Age = [options.Model.AgeMin, options.Model.AgeMax];
        MinMax_Model_Metalicity = [options.Model.MetallicityMin, options.Model.MetallicityMax];
        Data_Medium = options.Input.DataMediumAir;
        Flux_Scale_Factor = options.Input.FluxScaleFactor;
        Write_Results = options.Input.WriteResults;
        Milky_Way_Reddening = options.Dust.MilkyWayReddening;
        HPF_Mode = options.Dust.HpfMode.Equals("on", StringComparison.OrdinalIgnoreCase);
        Dust_Law = options.Dust.DustLaw.ToLowerInvariant() switch
        {
            "allen" => (byte)1,
            "prevot" => (byte)2,
            _ => (byte)0
        };
        Downgrade_models = options.Model.DowngradeModels;
        N_Masked_Amstrongs = options.Masking.MaskedAngstroms;
        Redshift = options.Input.Redshift;
        RA_DEC = [options.Input.Ra, options.Input.Dec];
        Velocity_Dispersion = options.Input.VelocityDispersion;
        Instrument_Resolution = options.Input.InstrumentResolution;
    }

    // SETUP
    public static void Setup(GPU gpu, FireflyOptions? options = null)
    {
        if (options is not null)
        {
            var models = Path.Combine(PathOfProgram, options.Paths.StellarPopulationModels.Replace('/', Path.DirectorySeparatorChar));
            var dust = Path.Combine(PathOfProgram, options.Paths.DustMaps.Replace('/', Path.DirectorySeparatorChar));
            Environment.SetEnvironmentVariable("FALCON_MODELS", models);
            Environment.SetEnvironmentVariable("FALCON_DUST", dust);
        }
        GetModelData(gpu, options);
        GetDustData(options);
        PreInitialiseDownGrade(gpu);
    }

    // SETUP (legacy)
    public static void Setup(GPU gpu) => Setup(gpu, null);

    // CONFIG PRE-INITIALISE DATA
    private static void GetModelData(GPU gpu, FireflyOptions? options)
    {
        string modelsDir = Environment.GetEnvironmentVariable("FALCON_MODELS")
            ?? Path.Combine(PathOfProgram, "assets", "stellar-population-models");
        string fileName = Path.Combine(modelsDir, "MaStar_SSP.h5");
        long id = Hdf5.OpenFile(fileName, true);

        // Load Data For Ma-Star Models
        if (Model_Key % 2 == 0)
        {
            try
            {
                // READ IN DATA
                FireflyConstants.r_model = (double[])((Hdf5.ReadDatasetToArray<double>(id, "r_model")).Result);

                FireflyConstants.t = (float[])((Hdf5.ReadDatasetToArray<float>(id, "t")).Result);
                FireflyConstants.Z = (float[])((Hdf5.ReadDatasetToArray<float>(id, "Z")).Result);
                FireflyConstants.s = (float[])((Hdf5.ReadDatasetToArray<float>(id, "s")).Result);

                FireflyConstants.wavelength = new Vector(gpu, (float[])((Hdf5.ReadDatasetToArray<float>(id, "wavelength")).Result));

                // READ IN FLUXGRID
                float[,,,] fluxgrid;
                if (Model_Key == 0b00001_010) // MaStar-Th
                {
                    fluxgrid = (float[,,,])((Hdf5.ReadDatasetToArray<float>(id, "fluxgrid_Th")).Result);
                }
                else if (Model_Key == 0b00010_010) // MaStar-E
                {
                    fluxgrid = (float[,,,])((Hdf5.ReadDatasetToArray<float>(id, "fluxgrid_E")).Result);
                }
                else
                {
                    throw new Exception("MaStar Model Flavour Error Please choose between MaStar-Th and MaStar-E");
                }

                Hdf5.CloseFile(id);


                // SET FLUXGRID 
                FireflyConstants.fluxgrid = fluxgrid;


                // SET SLOPE AND SIDX
                FireflyConstants.slope = IMF switch
                {
                    // kr
                    0 => 1.3f,
                    // ss
                    1 => 2.35f,
                    _ => throw new Exception("Unrecognised IMF"),
                };
                FireflyConstants.sidx = Array.IndexOf(FireflyConstants.s, FireflyConstants.slope);

            }
            catch (Exception)
            {
                throw new Exception("MaStar SSP Data not found");
            }

            return;
        }


        if (Model_Key % 2 != 0)
        {

            switch (Model_Key)
            {
                case 0b00001_001:
                    FireflyConstants.r_model = [2.55d];
                    break;
                case 0b00010_001:
                    FireflyConstants.r_model = [3.40d];
                    break;
                case 0b00100_001:
                    FireflyConstants.r_model = [0.55d];
                    break;
                case 0b01000_001:
                    FireflyConstants.r_model = [0.10d];
                    break;
            }

            return;
        }
    }

    /// <summary>
    /// Computes the Data Invariant section of DownGrade Function Outputting sres
    /// </summary>
    private static void PreInitialiseDownGrade(GPU gpu)
    {

        if (FireflyConstants.r_model.Length == 1)
        {

            FireflyConstants.sres = new Vector(gpu, new float[FireflyConstants.wavelength.Length]);

            for (int i = 0; i < FireflyConstants.r_model.Length; i++)
            {
                FireflyConstants.sres[i] = (float)(FireflyConstants.wavelength[i] * (1f / FireflyConstants.r_model[0]));
            }
            return;
        }

        FireflyConstants.sres = new(gpu, new float[FireflyConstants.r_model.Length]);
        FireflyConstants.sres.SyncCPU();
        for (int i = 0; i < FireflyConstants.r_model.Length; i++)
        {
            FireflyConstants.sres[i] = (float)FireflyConstants.r_model[i];
        }

        FireflyConstants.sres.UpdateCache();

        return;

    }

    private static void GetDustData(FireflyOptions? options)
    {
        string dustDir = Environment.GetEnvironmentVariable("FALCON_DUST")
            ?? Path.Combine(PathOfProgram, "assets", "dust-maps");
        string fileName = Path.Combine(dustDir, "dust.h5");
        long id = Hdf5.OpenFile(fileName, true);

        // Read Data
        FireflyConstants.ngp_dust = (float[,])((Hdf5.ReadDatasetToArray<float>(id, "ngp")).Result);
        FireflyConstants.sgp_dust = (float[,])((Hdf5.ReadDatasetToArray<float>(id, "sgp")).Result);
        return;
    }

}
