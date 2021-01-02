using System;
using System.Collections.Generic;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    public class Config
    {

        // CONFIG OF PHYSICS
        #region

        // Cosmology Values to use - DEFAULT = 0 : Planck15
        public static byte Cosmology { get; private set; } = 0;
        // Redshift Value
        public static float Redshift { get; private set; } = 0.021275453f;
        // Right Ascension and Declination
        public static float[] RA_DEC { get; private set; } = new float[2] { 145.89219f, 0.059372f };
        // Velocity Dispersion in km/s
        public static float Velocity_Dispersion { get; private set; } = 135.89957f;
        // Instrument Resolution !!errorNES!!
        public static float Instrument_Resolution { get; private set; } = 2000f;
        // Number masked amstrongs !!errorNES!!
        public static ushort N_Masked_Amstrongs { get; private set; } = 20;
        // Choose Model - DEFAULT = 0 : MaStar, 1 : m11
        public static byte Model_Key { get; private set; } = 0;
        // Choose Model Flavour 
        // FOR Model 0 - DEFAULT = 0 : E-MaStar, 1 : Th-MaStar
        // FOR Model 1 - DEFAULT = 0 : MILES, 1 : STELIB, 2 : ELODIE, 3 : MARCS (kr imf only)
        public static byte Model_Flavour { get; private set; } = 0;
        // Choose IMF : DEFAULT 0 : Kroupa, 1 : Salpeter
        public static byte IMF { get; private set; } = 0;
        // Min and Max Age of Models
        public static float[] MinMax_Model_Age { get; private set; } = new float[2] { 0, Constants.AoU };
        // Min and Max Metalicity of Models
        public static float[] MinMax_Model_Metalicity { get; private set; } = new float[2] { -3f, 3f };
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

        #endregion

        // CONFIG OF AI

    }
}
