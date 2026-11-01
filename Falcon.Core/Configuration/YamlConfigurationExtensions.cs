using Falcon.Core.Configuration;
using YamlDotNet.Serialization;
using YamlDotNet.Serialization.NamingConventions;

namespace Falcon.Core.Configuration;

public static class YamlConfigurationExtensions
{
    public static FireflyOptions LoadFireflyYaml(string path)
    {
        using var reader = File.OpenText(path);
        return Deserialize(reader);
    }

    public static FireflyOptions LoadFireflyYamlFromString(string yaml) =>
        Deserialize(new StringReader(yaml));

    static FireflyOptions Deserialize(TextReader reader)
    {
        var deserializer = new DeserializerBuilder()
            .WithNamingConvention(UnderscoredNamingConvention.Instance)
            .IgnoreUnmatchedProperties()
            .Build();
        return deserializer.Deserialize<FireflyOptions>(reader)
            ?? throw new InvalidOperationException("Failed to deserialize YAML.");
    }

    public static void ApplyToLegacy(FireflyOptions options, string basePath) =>
        LegacyConfig.ApplyFrom(options, basePath);
}
