namespace Falcon.Core.Services.Sfr;

public interface IStarFormationRateService
{
    double EstimateSfr(double[] agePdf, double[] ageBins);
}

public sealed class StarFormationRateService : IStarFormationRateService
{
    public double EstimateSfr(double[] agePdf, double[] ageBins)
    {
        if (agePdf.Length != ageBins.Length) throw new ArgumentException("PDF and age bins length mismatch.");
        double sfr = 0;
        for (int i = 0; i < agePdf.Length; i++)
            sfr += agePdf[i] * (ageBins[i] < 0.1 ? 1 : 0);
        return sfr;
    }
}
