using System;
using System.Collections.Generic;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    
    public struct Constants
    {
        // Example
        // readonly public static float pi = 3.14f;
        
        // FORMAT : { Emission Name : Emission Values }
        readonly public static Dictionary<string, float[]> Emission_lines = new Dictionary<string, float[]>
        {
            { "He-II" , new float[2] { 3202.15f, 4685.74f } },
            { "Ne-V" , new float[2] { 3345.81f, 3425.81f } },
            { "O-II" , new float[2] { 3726.03f, 3728.73f } },
            { "Ne-III" , new float[2] { 3868.69f, 3967.40f } },
            { "H-ζ" , new float[1] { 3889.05f } },
            { "H-ε" , new float[1] { 3970.07f } },
            { "H-δ" , new float[1] { 4101.73f } },
            { "H-γ" , new float[1] { 4340.46f } },
            { "O-III" , new float[3] { 4363.15f, 4958.83f, 5006.77f } },           
            { "Ar-IV" , new float[2] { 4711.30f, 4740.10f } },
            { "H-β" , new float[1] { 4861.32f } },
            { "N-I" , new float[2] { 5197.90f, 5200.39f } },
            { "He-I" , new float[1] { 5875.60f } },
            { "O-I" , new float[2] { 6300.20f, 6363.67f } },
            { "N-II" , new float[2] { 6547.96f, 6583.34f } },
            { "H-α" , new float[1] { 6562.80f } },
            { "S-II" , new float[2] { 6716.31f, 6730.68f } },
            { "Ar-III" , new float[1] { 7135.67f } },
        };


        // AGE OF UNIVERSE in YEARS
        readonly public static float AoU = 14.7e9f; // !!errorVaL!!
        // The speed of light in meters per second
        readonly public static float c = 299792458f;



        // Conversion Factor From Degrees to Radians
        readonly public static float Deg2RadFactor = (float)Math.PI / 180f;
        // Conversion Factor From Degrees to Radians
        readonly public static float Rad2DegFactor = 180f / (float)Math.PI;
        // Conversion Factor From Degrees to Radians
        readonly public static float TwoPi = (float)Math.PI * 2f;



        // RA(radians),Dec(radians),distance(kpc) of Galactic center in J2000
        readonly public static Vector Galactic_Center_Equatorial = 
            new Vector(
                new float[3] { 
                    UtilityMethods.Degree2Radians(266.40510f), 
                    UtilityMethods.Degree2Radians(-28.936175f), 
                    8.33f }, 1
                );
        // RA(radians),Dec(radians) of Galactic Northpole in J2000
        readonly public static Vector Galactic_Northpole_Equatorial =
            new Vector(
                new float[2] {
                    UtilityMethods.Degree2Radians(192.859508f),
                    UtilityMethods.Degree2Radians(27.128336f)
                    }, 1
                );

        // Cosmology Parameters
        readonly public static double H0units_to_invs = 3.240779289469756e-20f;
        readonly public static double sec_to_Gyr = 3.168808781402895e-17f;
        readonly public static double critdens_const = 1788445.339869672f; // g/cm^3
        readonly public static double arcsec_in_radians = 4.84813681109536e-06f;
        readonly public static double arcmin_in_radians = 0.0002908882086657216f;
        readonly public static double a_B_c2 = 8.418013525010775e-36f;
        readonly public static double kB_evK = 8.617333262145179e-05f; // eV/K

    }
}
