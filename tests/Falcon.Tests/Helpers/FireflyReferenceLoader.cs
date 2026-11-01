using System.Text.Json;

namespace Falcon.Tests.Helpers;

public static class FireflyReferenceLoader
{
    public static JsonDocument LoadGolden(string fileName)
    {
        string path = Path.Combine(AppContext.BaseDirectory, "TestData", fileName);
        return JsonDocument.Parse(File.ReadAllText(path));
    }

    public static double[] GetDoubleArray(JsonDocument doc, string property)
    {
        var arr = doc.RootElement.GetProperty(property);
        return arr.EnumerateArray().Select(e => e.GetDouble()).ToArray();
    }

    public static double GetDouble(JsonDocument doc, string property) =>
        doc.RootElement.GetProperty(property).GetDouble();
}
