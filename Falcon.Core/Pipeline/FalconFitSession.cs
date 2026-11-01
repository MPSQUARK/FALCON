namespace Falcon.Core.Pipeline;

/// <summary>
/// Coordinates cancellation for in-flight spectral fits when the host shuts down.
/// Batch cancellation is owned by the caller (e.g. <c>RunBatchAsync</c> linked token).
/// </summary>
public sealed class FalconFitSession : IDisposable
{
    readonly CancellationTokenSource _shutdownCts = new();

    public bool IsShutdown => _shutdownCts.IsCancellationRequested;

    public CancellationToken ShutdownToken => _shutdownCts.Token;

    public CancellationToken GetFitToken(CancellationToken callerToken = default)
    {
        if (!callerToken.CanBeCanceled)
            return _shutdownCts.Token;

        return CancellationTokenSource.CreateLinkedTokenSource(_shutdownCts.Token, callerToken).Token;
    }

    public void RequestShutdown()
    {
        try { _shutdownCts.Cancel(); }
        catch (ObjectDisposedException) { }
    }

    public void Dispose()
    {
        RequestShutdown();
        _shutdownCts.Dispose();
    }
}
