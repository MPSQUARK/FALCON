using BAVCL;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Dust;
using FluentAssertions;

namespace Falcon.Tests.Unit;

public class SfdDustTests : IDisposable
{
    private readonly GPU _gpu = FalconGpuHost.Create();
    private readonly DustService _dust;

    public SfdDustTests() => _dust = new DustService(_gpu);

    public void Dispose() => _gpu.Dispose();

    private static float[,] BuildGrid(float v00, float v01, float v10, float v11) =>
        new[,] { { v00, v01 }, { v10, v11 } };

    [Fact]
    public void GetSfdEbv_NorthernHemisphere_UsesNgpGrid()
    {
        var ngp = BuildGrid(0.01f, 0.02f, 0.03f, 0.04f);
        var sgp = BuildGrid(0.10f, 0.20f, 0.30f, 0.40f);
        var ebv = _dust.GetSfdEbv(longitudeDeg: 0, latitudeDeg: 45, ngp, sgp);
        ebv.Should().BeApproximately(0.01, 0.05);
    }

    [Fact]
    public void GetSfdEbv_SouthernHemisphere_UsesSgpGrid()
    {
        var ngp = BuildGrid(0.01f, 0.02f, 0.03f, 0.04f);
        var sgp = BuildGrid(0.10f, 0.20f, 0.30f, 0.40f);
        var ebv = _dust.GetSfdEbv(longitudeDeg: 0, latitudeDeg: -45, ngp, sgp);
        ebv.Should().BeApproximately(0.3, 0.05);
    }

    [Theory]
    [InlineData(0, 30, 0)]
    [InlineData(90, 30, 0)]
    [InlineData(180, 30, 1)]
    [InlineData(270, 30, 1)]
    public void GetSfdEbv_LongitudeWrapsOntoGrid(double lon, double lat, int expectedCol)
    {
        var ngp = BuildGrid(1f, 2f, 3f, 4f);
        var sgp = BuildGrid(5f, 6f, 7f, 8f);
        var ebv = _dust.GetSfdEbv(lon, lat, ngp, sgp);
        var expected = expectedCol == 0 ? 1f : 2f;
        ebv.Should().BeApproximately(expected, 1.0);
    }

    [Fact]
    public void GetSfdEbv_Equator_SelectsBySignOfLatitude()
    {
        var ngp = BuildGrid(0.05f, 0.06f, 0.07f, 0.08f);
        var sgp = BuildGrid(0.15f, 0.16f, 0.17f, 0.18f);
        var north = _dust.GetSfdEbv(45, 0.01, ngp, sgp);
        var south = _dust.GetSfdEbv(45, -0.01, ngp, sgp);
        north.Should().NotBe(south);
    }

    [Theory]
    [InlineData(45)]
    [InlineData(89)]
    public void GetSfdEbv_NorthernLatitude_ClampsRowIndex(double lat)
    {
        var ngp = BuildGrid(0.1f, 0.2f, 0.3f, 0.4f);
        var sgp = BuildGrid(0.5f, 0.6f, 0.7f, 0.8f);
        var ebv = _dust.GetSfdEbv(0, lat, ngp, sgp);
        ebv.Should().BeInRange(0.0, 0.5);
    }
}
