using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    public class Spectra
    {
        public Spectra(string path,
                       bool milky_Way_Reddening = true,
                       bool hPF_Mode = true,
                       ushort n_Masked_Amstrongs = 20)
        {
            this.Path = path;
            this.Milky_Way_Reddening = milky_Way_Reddening;
            this.HPF_Mode = hPF_Mode;
            this.N_Masked_Amstrongs = n_Masked_Amstrongs;
        }

        // constructor variables
        public string Path { get; private set; }
        public bool Milky_Way_Reddening { get; private set; }
        public bool HPF_Mode { get; private set; }
        public ushort N_Masked_Amstrongs { get; private set; }

        // Initialisation variables
        public Vector Wavelength { get; private set; }
        public Vector Flux { get; private set; }
        public Vector Error { get; private set; }
        public float Redshift { get; private set; }
        public float[] RA_DEC { get; private set; }
        public float Velocity_Dispersion { get; private set; }
        public float Instrument_Resolution { get; private set; }

        public float Distance_Luminosity { get; private set; }
        public Vector Bad_Flags { get; private set; }
        public Vector Restframe_Wavelength { get; private set; }
        // !!errorNES!!
        public byte Trust_Flag { get; private set; }
        // !!errorNES!!
        public byte ObjID { get; private set; }
        // Bad Data ?? 
        // mask emission lines
        public float ebv_MW { get; private set; }




        public static Spectra InitialiseSpectraParameters(AcceleratorId acceleratorId, Spectra spectrum, Vector Data, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
        {
            spectrum.Wavelength = Vector.AccessSlice(acceleratorId, Data, 0, 'c');
            spectrum.Flux = Vector.AccessSlice(acceleratorId, Data, 1, 'c');
            spectrum.Error = Vector.AccessSlice(acceleratorId, Data, 2, 'c');
            spectrum.Redshift = Redshift;
            spectrum.RA_DEC = RA_DEC;
            spectrum.Velocity_Dispersion = Velocity_Disp;
            spectrum.Instrument_Resolution = instrument_resolution;

            // CALCULATE LUMINOSITY DISTANCE in CM
            spectrum.Bad_Flags = Vector.Fill(acceleratorId, 1, spectrum.Wavelength.Value.Length);
            spectrum.Restframe_Wavelength = Vector.ScalarOperation(acceleratorId, spectrum.Wavelength, (1 + spectrum.Redshift), '/');
            spectrum.Trust_Flag = 1;
            spectrum.ObjID = 0;

            // Remove Bad data from the spectrum
            spectrum.MaskEmissionlines();

            if (spectrum.Milky_Way_Reddening)
            {
                spectrum.ebv_MW = spectrum.GetDustRADEC(spectrum.RA_DEC, "ebv");
            }
            else
            {
                spectrum.ebv_MW = 0f;
            }


            return spectrum;
        }

        private void MaskEmissionlines()
        {
            
        }

        private float GetDustRADEC(float[] RADEC, string dustmap, bool interpolate = true)
        {

            float[] l_b = Eq2Gal(
                ProcessDataMethods.Degree2Radians(RADEC[0]), 
                ProcessDataMethods.Degree2Radians(RADEC[1])
                );

            return Get_SFD_dust(l_b[0], l_b[1]);
        }

        private float[] Eq2Gal(float ra, float dec)
        {

            float Cos_dec_Mult_Cos_GNEV0 = (float)(Math.Cos(dec) * Math.Cos(ra - Constants.Galactic_Northpole_Equatorial.Value[0]));

            float b = (float)Math.Asin(
                Math.Sin(dec) * Math.Sin(Constants.Galactic_Northpole_Equatorial.Value[1]) + 
                Math.Cos(Constants.Galactic_Northpole_Equatorial.Value[1]) * Cos_dec_Mult_Cos_GNEV0
                );

            float l = (float)Math.Atan2(
                Math.Sin(dec) * Math.Cos(Constants.Galactic_Northpole_Equatorial.Value[1]) - 
                Math.Sin(Constants.Galactic_Northpole_Equatorial.Value[1]) * Cos_dec_Mult_Cos_GNEV0, 
                Math.Cos(dec) * Math.Sin(ra - Constants.Galactic_Northpole_Equatorial.Value[0])
                ) + ProcessDataMethods.Degree2Radians(33f);

            if (l < 0)
            {
                l += Constants.TwoPi;
            }

            l %= Constants.TwoPi;


            return new float[2] { l * Constants.Rad2DegFactor, b * Constants.Rad2DegFactor };
        }

        private float Get_SFD_dust(float logitude, float latitude, string dustmap = "ebv", bool interpolate = true)
        {
            // get path to dust map



            return 0f;
        }

    }
}
