namespace Falcon.Core.Services.MassLoss;

/// <summary>Inverse Squared Distance interpolation from firefly_estimations_3d.py</summary>
public sealed class MassLossEstimationService
{
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[] _v;

    public MassLossEstimationService(double[] dataX, double[] dataY, double[] dataZ)
    {
        _x = dataX;
        _y = dataY;
        _v = dataZ;
    }

    public double Estimate(double x, double y) => Isd(x, y);

    private double Isd(double x, double y)
    {
        var d = _x.Zip(_y, (xi, yi) => Math.Sqrt((x - xi) * (x - xi) + (y - yi) * (y - yi))).ToArray();
        double dMin = d.Min();
        if (dMin > 0)
        {
            var weights = d.Select(di => 1.0 / (di * di)).ToArray();
            double wSum = weights.Sum();
            return _v.Zip(weights, (vi, wi) => vi * wi / wSum).Sum();
        }
        int idx = Array.IndexOf(d, dMin);
        return _v[idx];
    }
}

public sealed class MassLossResult
{
    public required double LivingStarsMass { get; init; }
    public required double RemnantMass { get; init; }
    public required double WhiteDwarfMass { get; init; }
    public required double NeutronStarMass { get; init; }
    public required double BlackHoleMass { get; init; }
    public required double EjectaMass { get; init; }
    public required double StellarMass { get; init; }
}
