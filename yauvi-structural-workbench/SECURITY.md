# Security policy

The workbench is local-only. It binds to loopback, checks Host and Origin,
requires a session token for mutations, invokes registered commands only, and
does not expose a filesystem browser, arbitrary URL fetcher, or external upload.

Report path traversal, command injection, arbitrary-file access, unsafe archive
expansion, source-host bypass, checksum substitution, and silent scientific
corruption privately once the public repository's security channel exists.

External tools and provider services retain their own security policies.
