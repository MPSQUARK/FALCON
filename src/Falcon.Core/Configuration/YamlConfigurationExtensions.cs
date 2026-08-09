using Falcon.Core.Configuration;
using YamlDotNet.Serialization;
using YamlDotNet.Serialization.NamingConventions;

namespace Falcon.Core.Configuration;

public static class YamlConfigurationExtensions
{
    public static FireflyOptions LoadFireflyYaml(string path)
    {
        using var reader = File.OpenText(path);
        var deserializer = new DeserializerBuilder()
            .WithNamingConvention(UnderscoredNamingConvention.Instance)
            .IgnoreUnmatchedProperties()
            .Build();
        return deserializer.Deserialize<FireflyOptions>(reader)
            ?? throw new InvalidOperationException($"Failed to deserialize YAML: {path}");
    }

    public static void ApplyToLegacy(FireflyOptions options, string basePath) =>
        LegacyConfig.ApplyFrom(options, basePath);
}
