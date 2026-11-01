using System.Globalization;
using Serilog.Core;
using Serilog.Events;
using Serilog.Formatting.Display;
using Serilog.Formatting;

namespace Falcon.UI.Services;

/// <summary>Serilog sink exposing the latest fit status line for the Run page.</summary>
public sealed class RunStatusSink : ILogEventSink
{
    static readonly MessageTemplateTextFormatter Formatter =
        new("{Timestamp:HH:mm:ss} [{Level:u3}] {Message:lj}{NewLine}{Exception}", CultureInfo.InvariantCulture);

    string _lastMessage = "";

    public event Action? Changed;

    public string LastMessage => _lastMessage;

    public void Emit(LogEvent logEvent)
    {
        using var writer = new StringWriter();
        Formatter.Format(logEvent, writer);
        _lastMessage = writer.ToString().TrimEnd();
        Changed?.Invoke();
    }
}
