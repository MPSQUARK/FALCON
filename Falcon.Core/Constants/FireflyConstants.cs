using BAVCL;

namespace Falcon.Core.Constants;


public static class FireflyConstants
{
    /* READ ONLY
     */

    // SPECTROSCOPY
    //
    public static readonly Dictionary<string, float[]> Emission_lines = new Dictionary<string, float[]>
    {
        { "He-II" , new float[2] { 3202.15f, 4685.74f } },
        { "Ne-V" , new float[2] { 3345.81f, 3425.81f } },
        { "O-II" , new float[2] { 3726.03f, 3728.73f } },
        { "Ne-III" , new float[2] { 3868.69f, 3967.40f } },
        { "H-zeta" , new float[1] { 3889.05f } },
        { "H-e" , new float[1] { 3970.07f } },
        { "H-d" , new float[1] { 4101.73f } },
        { "H-gamma" , new float[1] { 4340.46f } },
        { "O-III" , new float[3] { 4363.15f, 4958.83f, 5006.77f } },           
        { "Ar-IV" , new float[2] { 4711.30f, 4740.10f } },
        { "H-�" , new float[1] { 4861.32f } },
        { "N-I" , new float[2] { 5197.90f, 5200.39f } },
        { "He-I" , new float[1] { 5875.60f } },
        { "O-I" , new float[2] { 6300.20f, 6363.67f } },
        { "N-II" , new float[2] { 6547.96f, 6583.34f } },
        { "H-a" , new float[1] { 6562.80f } },
        { "S-II" , new float[2] { 6716.31f, 6730.68f } },
        { "Ar-III" , new float[1] { 7135.67f } },
    };



    // UNIVERSAL CONSTANTS

    // AGE OF UNIVERSE in YEARS
    public static readonly float AoU = 14.7e9f; // !!errorVaL!!
    // The speed of light in meters per second
    public static readonly float c = 299792458f;
    public static readonly float c_kms = 299792.458f;
    public static readonly float sig2FWHM = 2f * MathF.Sqrt(2f * MathF.Log(2));
    public static readonly float c_div_sig2 = c_kms / sig2FWHM;
    // length of 1 sidereal year in seconds
    public static readonly float sidereal_year = 31558175.779f;


    // CONVERSION FACTORS

    // Conversion Factor From Degrees to Radians
    public static readonly float Deg2RadFactor = MathF.PI / 180f;
    // Conversion Factor From Degrees to Radians
    public static readonly float Rad2DegFactor = 180f / MathF.PI;
    // Conversion Factor From Degrees to Radians
    public static readonly float TwoPi = MathF.PI * 2f;
    // Conversion Factor From Radians to arcseconds
    public static readonly float rad2arcsec = 3600f * Rad2DegFactor;


    // OBSERVATIONAL ASTRONOMY PARAMETERS

    // RA(radians),Dec(radians),distance(kpc) of Galactic center in J2000
    public static readonly float[] Galactic_Center_Equatorial = 
            new float[3] { 
                266.40510f * Deg2RadFactor, 
                -28.936175f * Deg2RadFactor, 
                8.33f 
            };
    // RA(radians),Dec(radians) of Galactic Northpole in J2000
    public static readonly float[] Galactic_Northpole_Equatorial =
            new float[2] {
                192.859508f * Deg2RadFactor,
                27.128336f * Deg2RadFactor
            };



    // Cosmology Parameters 
    //
    public static readonly double H0units_to_invs = 3.240779289469756e-20f;
    public static readonly double sec_to_Gyr = 3.168808781402895e-17f;
    public static readonly double critdens_const = 1788445.339869672f; // g/cm^3
    public static readonly double arcsec_in_radians = 4.84813681109536e-06f;
    public static readonly double arcmin_in_radians = 0.0002908882086657216f;
    public static readonly double a_B_c2 = 8.418013525010775e-36f;
    public static readonly double kB_evK = 8.617333262145179e-05f; // eV/K


    /* GLOBAL DATA
     */

    // Models - DATA

    // MaStar
    public static double[] r_model { get; set; }
    public static float[] t { get; set; }
    public static float[] Z { get; set; }
    public static float[] s { get; set; }
    public static Vector wavelength { get; set; }
    /// <summary>CPU copy of MaStar wavelength grid (GPU vector may be invalidated during downgrade).</summary>
    public static float[] wavelengthCpu { get; set; } = [];
    public static float[,,,] fluxgrid { get; set; } // Line 326-329 - firefly models



    // Models - Other
    //
    public static float slope { get; set; } = 0f;
    public static int sidx { get; set; } = 0;
    public static Vector sres { get; set; }

    // Dust Maps
    public static float[,] ngp_dust { get; set; }
    public static float[,] sgp_dust { get; set; }

}
