using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MachineLearningSpectralFittingCode
{
    public class SPModel
    {
        int
            velocity_dispersion_r { get; set; }
        int
            fit_per_iteration_cap {get; set;}
        
        List<double> 
            delta_lamdba_lib { get; set; }

        public void InitialiseSPModel(float velocity_disp)
        {
            this.delta_lamdba_lib = new List<double>();
            switch (Program.config.Model_Key)
            {
                case 0b00001_001:
                    this.delta_lamdba_lib.Add(2.55d);
                    break;
                case 0b00010_001:
                    this.delta_lamdba_lib.Add(3.40d);
                    break;
                case 0b00100_001:
                    this.delta_lamdba_lib.Add(0.55d);
                    break;
                case 0b01000_001:
                    this.delta_lamdba_lib.Add(0.10d);
                    break;


                case 0b00001_010:
                    this.delta_lamdba_lib.AddRange(Constants.MaStarSSP);
                    break;
                case 0b00010_010:                    
                    this.delta_lamdba_lib.AddRange(Constants.MaStarSSP);
                    break;

                default:
                    throw new Exception("Incorrect Model key");
            }
            this.velocity_dispersion_r = (int)(MathF.Round(velocity_disp / 5f) * 5f);
            
        }



    }
}
