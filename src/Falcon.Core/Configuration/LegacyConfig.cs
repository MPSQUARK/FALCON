using BAVCL;
using BAVCL.Core;
using Falcon.Core.Constants;
using Falcon.Core.Services.Cosmology;
using Falcon.Core.Services.IO;

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
    // Set Parameters for dust determination: on | hpf_only | off
    public static string HpfMode { get; private set; } = "on";
    public static bool HPF_Mode => HpfMode.Equals("on", StringComparison.OrdinalIgnoreCase);
    public static bool HpfOnlyMode => HpfMode.Equals("hpf_only", StringComparison.OrdinalIgnoreCase);
    // Set Dust Law - DEFAULT 0 : Calzetti, 1 : Allen, 2 : prevot
    public static byte Dust_Law { get; private set; } = 0;
    public static bool Downgrade_models { get; private set; } = true;

    public static string FittingEngine { get; private set; } = "firefly";
    public static bool LegacyGpuFallback { get; private set; } = true;
    public static int MaxIterations { get; private set; } = 10;
    public static int FitPerIterationCap { get; private set; } = 1000;
    public static int PdfSampling { get; private set; } = 200;
    public static float MaxEbv { get; private set; } = 1.0f;
    public static int NumDustVals { get; private set; } = 200;
    public static float DustSmoothingLength { get; private set; } = 100f;
    public static int TopResultsCount { get; private set; } = 3;
    public static int MaxSspComponents { get; private set; } = 5;
    public static bool UseGpuEvaluator { get; private set; } = false;
    public static bool ExperimentalMixSearch { get; private set; } = false;
    public static int ExperimentalMixSearchMaxEvals { get; private set; } = 64;
    public static int MaxParallelFits { get; private set; }
    public static double[] WaveLimits { get; private set; } = [0, 99_999_990];
    public static float[] EmissionLines { get; private set; } = [6562.80f, 6583.34f, 6547.96f, 4861.32f, 5006.77f, 4958.83f];

    #endregion

    public static void ApplyFrom(FireflyOptions options, string basePath)
    {
        PathOfProgram = basePath.EndsWith(Path.DirectorySeparatorChar) ? basePath : basePath + Path.DirectorySeparatorChar;
        Model_Key = options.Model.ModelKey;
        IMF = options.Model.Imf.Equals("Salpeter", StringComparison.OrdinalIgnoreCase) ? (byte)1 : (byte)0;
        MinMax_Model_Age = [NormalizeAgeGyr(options.Model.AgeMin), NormalizeAgeGyr(options.Model.AgeMax)];
        MinMax_Model_Metalicity = [options.Model.MetallicityMin, options.Model.MetallicityMax];
        Data_Medium = options.Input.DataMediumAir;
        Flux_Scale_Factor = options.Input.FluxScaleFactor;
        Write_Results = options.Input.WriteResults;
        Milky_Way_Reddening = options.Dust.MilkyWayReddening;
        HpfMode = options.Dust.HpfMode.ToLowerInvariant();
        Dust_Law = options.Dust.DustLaw.ToLowerInvariant() switch
        {
            "allen" => (byte)1,
            "prevot" => (byte)2,
            _ => (byte)0
        };
        Downgrade_models = options.Model.DowngradeModels;
        N_Masked_Amstrongs = options.Masking.MaskedAngstroms;
        EmissionLines = options.Masking.EmissionLines;
        FittingEngine = options.Fitting.Engine.ToLowerInvariant();
        LegacyGpuFallback = options.Fitting.LegacyGpuFallback;
        MaxIterations = options.Fitting.MaxIterations;
        FitPerIterationCap = options.Fitting.FitPerIterationCap;
        PdfSampling = options.Fitting.PdfSampling;
        MaxEbv = options.Fitting.MaxEbv;
        NumDustVals = options.Fitting.NumDustVals;
        DustSmoothingLength = options.Fitting.DustSmoothingLength;
        TopResultsCount = Math.Max(1, options.Fitting.TopResultsCount);
        MaxSspComponents = Math.Max(1, options.Fitting.MaxSspComponents);
        UseGpuEvaluator = options.Fitting.UseGpuEvaluator;
        ExperimentalMixSearch = options.Fitting.ExperimentalMixSearch;
        ExperimentalMixSearchMaxEvals = Math.Max(8, options.Fitting.ExperimentalMixSearchMaxEvals);
        MaxParallelFits = Math.Max(0, options.Fitting.MaxParallelFits);
        WaveLimits = options.Fitting.WaveLimits.Length >= 2
            ? [options.Fitting.WaveLimits[0], options.Fitting.WaveLimits[1]]
            : [0, 99_999_990];
        Redshift = options.Input.Redshift;
        RA_DEC = [options.Input.Ra, options.Input.Dec];
        Velocity_Dispersion = options.Input.VelocityDispersion;
        Instrument_Resolution = options.Input.InstrumentResolution;
    }

    // SETUP
    private static byte _loadedModelKey;
    private static byte _loadedImf;
    private static string? _loadedModelsPath;
    private static string? _loadedDustPath;
    private static bool _modelDataLoaded;
    private static bool _dustDataLoaded;

    public static void Setup(GPU gpu, FireflyOptions? options = null)
    {
        string modelsDir = Environment.GetEnvironmentVariable("FALCON_MODELS")
            ?? Path.Combine(PathOfProgram, "assets", "stellar-population-models");
        string dustDir = Environment.GetEnvironmentVariable("FALCON_DUST")
            ?? Path.Combine(PathOfProgram, "assets", "dust-maps");
        if (options is not null)
        {
            modelsDir = Path.Combine(PathOfProgram, options.Paths.StellarPopulationModels.Replace('/', Path.DirectorySeparatorChar));
            dustDir = Path.Combine(PathOfProgram, options.Paths.DustMaps.Replace('/', Path.DirectorySeparatorChar));
            Environment.SetEnvironmentVariable("FALCON_MODELS", modelsDir);
            Environment.SetEnvironmentVariable("FALCON_DUST", dustDir);
        }
        modelsDir = Path.GetFullPath(modelsDir);
        dustDir = Path.GetFullPath(dustDir);

        if (!_modelDataLoaded || _loadedModelKey != Model_Key || _loadedImf != IMF
            || !string.Equals(_loadedModelsPath, modelsDir, StringComparison.OrdinalIgnoreCase)
            || FireflyConstants.fluxgrid is null)
        {
            GetModelData(gpu, modelsDir);
            _loadedModelKey = Model_Key;
            _loadedImf = IMF;
            _loadedModelsPath = modelsDir;
            _modelDataLoaded = true;
        }

        if (!_dustDataLoaded
            || !string.Equals(_loadedDustPath, dustDir, StringComparison.OrdinalIgnoreCase)
            || FireflyConstants.ngp_dust is null)
        {
            GetDustData(dustDir);
            _loadedDustPath = dustDir;
            _dustDataLoaded = true;
        }

        PreInitialiseDownGrade(gpu);
    }

    // SETUP (legacy)
    public static void Setup(GPU gpu) => Setup(gpu, null);

    // CONFIG PRE-INITIALISE DATA
    private static void GetModelData(GPU gpu, string modelsDir)
    {
        string fileName = Path.Combine(modelsDir, "MaStar_SSP_v0.2.fits.gz");

        // Load Data For Ma-Star Models
        if (Model_Key % 2 == 0)
        {
            try
            {
                const int ages = 24;
                const int metallicities = 9;
                const int imfs = 8;
                const int wavelengths = 4563;
                int fluxHdu = Model_Key switch
                {
                    0b00001_010 => 3,
                    0b00010_010 => 4,
                    _ => throw new InvalidOperationException(
                        "MaStar Model Flavour Error: choose between MaStar-Th and MaStar-E.")
                };
                var modelImages = FitsAssetReader.ReadImages(fileName, 1, 2, fluxHdu);
                double[] parameters = modelImages[1];
                double[] waveAndResolution = modelImages[2];

                FireflyConstants.t = ReadParameterAxis(parameters, ages, metallicities, imfs, 0, 0);
                FireflyConstants.Z = ReadParameterAxis(parameters, ages, metallicities, imfs, 1, 1);
                FireflyConstants.s = ReadParameterAxis(parameters, ages, metallicities, imfs, 2, 2);
                FireflyConstants.r_model = waveAndResolution
                    .Skip(wavelengths)
                    .Take(wavelengths)
                    .ToArray();
                float[] wavelength = waveAndResolution
                    .Take(wavelengths)
                    .Select(value => (float)value)
                    .ToArray();
                FireflyConstants.wavelength = new Vector(gpu, wavelength);
                FireflyConstants.wavelengthCpu = FireflyConstants.wavelength.RetrieveReadOnlySpan().ToArray();

                double[] fluxValues = modelImages[fluxHdu];
                FireflyConstants.fluxgrid = ToFluxGrid(
                    fluxValues, ages, metallicities, imfs, wavelengths);

                // SET SLOPE AND SIDX
                FireflyConstants.slope = IMF switch
                {
                    0 => 1.3f,
                    1 => 2.35f,
                    _ => throw new Exception("Unrecognised IMF"),
                };
                FireflyConstants.sidx = Array.IndexOf(FireflyConstants.s, FireflyConstants.slope);

            }
            catch (Exception ex) when (ex is IOException or InvalidDataException or IndexOutOfRangeException)
            {
                throw new InvalidOperationException(
                    $"MaStar SSP FITS data not found or invalid: {fileName}", ex);
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
            FireflyConstants.sres = FireflyConstants.wavelength * (float)(1.0 / FireflyConstants.r_model[0]);
            return;
        }

        float[] sresValues = new float[FireflyConstants.r_model.Length];
        for (int i = 0; i < sresValues.Length; i++)
            sresValues[i] = (float)FireflyConstants.r_model[i];
        FireflyConstants.sres = new Vector(gpu, sresValues);

        return;

    }

    private static void GetDustData(string dustDir)
    {
        string ngpPath = Path.Combine(dustDir, "SFD_dust_4096_ngp.fits");
        string sgpPath = Path.Combine(dustDir, "SFD_dust_4096_sgp.fits");
        FireflyConstants.ngp_dust = ToDustMap(FitsAssetReader.ReadImage(ngpPath));
        FireflyConstants.sgp_dust = ToDustMap(FitsAssetReader.ReadImage(sgpPath));
    }

    private static float[] ReadParameterAxis(
        double[] values,
        int ages,
        int metallicities,
        int imfs,
        int component,
        int axis)
    {
        int length = axis switch
        {
            0 => ages,
            1 => metallicities,
            2 => imfs,
            _ => throw new ArgumentOutOfRangeException(nameof(axis))
        };
        var result = new float[length];
        for (int index = 0; index < length; index++)
        {
            int age = axis == 0 ? index : 0;
            int metallicity = axis == 1 ? index : 0;
            int imf = axis == 2 ? index : 0;
            int offset = (((age * metallicities) + metallicity) * imfs + imf) * 3 + component;
            result[index] = (float)values[offset];
        }
        return result;
    }

    private static float[,,,] ToFluxGrid(
        double[] values,
        int ages,
        int metallicities,
        int imfs,
        int wavelengths)
    {
        int expected = ages * metallicities * imfs * wavelengths;
        if (values.Length != expected)
            throw new InvalidDataException(
                $"MaStar flux image has {values.Length} values; expected {expected}.");

        var grid = new float[ages, metallicities, imfs, wavelengths];
        int offset = 0;
        for (int age = 0; age < ages; age++)
        for (int metallicity = 0; metallicity < metallicities; metallicity++)
        for (int imf = 0; imf < imfs; imf++)
        for (int wavelength = 0; wavelength < wavelengths; wavelength++)
            grid[age, metallicity, imf, wavelength] = (float)values[offset++];
        return grid;
    }

    private static float[,] ToDustMap(double[] values)
    {
        const int size = 4096;
        if (values.Length != size * size)
            throw new InvalidDataException(
                $"SFD dust image has {values.Length} values; expected {size * size}.");

        var map = new float[size, size];
        int offset = 0;
        for (int row = 0; row < size; row++)
        for (int column = 0; column < size; column++)
            map[row, column] = (float)values[offset++];
        return map;
    }

    /// <summary>MaStar ages are in Gyr (firefly convention). Values &gt;100 are treated as years.</summary>
    internal static float NormalizeAgeGyr(float age) =>
        age > 100f ? age / 1e9f : age;

}
