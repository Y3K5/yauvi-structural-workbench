# Opening the current Structural Workbench

The product is **YAUVI Structural Workbench**. Its browser is a local application,
with the page title and footer identifying the version and application build.
The default port is **8947**, keeping it distinct from the older broader YAUVI
portal on 8931. An explicit `--port` continues to work.

Installed distribution:

```sh
yauvi --workspace ./my-analysis workbench open
```

This opens the default browser and serves the local application. Keep the server
terminal open. Repeating the command opens the existing server only if its
application build, browser assets, analysis library and retrieval settings match.
An older or different server is reported, not terminated or silently reused.
Use `--label "Public preview" --port 8948` to distinguish another local copy.
The label appears in the page header, browser tab and About panel. Give each
copy its own `--workspace` so experimental inputs never enter a demo library.
`workbench serve` starts the same application without opening a browser.

For a development checkout, use a Python environment with the declared project
dependencies installed:

```sh
.venv/bin/python tools/launch_workbench.py --workspace ./my-analysis workbench open
```

This resolves application and engine source from that checkout, independent of
the caller's directory. It does not install dependencies or download updates.
Use `--allow-reference-fetch` only when enabling deliberate public-source
retrieval. Opening the portal does not retrieve protein files.

## Reading the identity

Click the version in the footer to open **About this window**. The semantic
version names the release; the content-based build distinguishes edits within
the same development version. The interface identity includes the local viewer,
membrane renderer, HTML, CSS and JavaScript. The application build also identifies
the packaged application and engine files. Absolute installation and analysis
paths are not returned by the build endpoint.

The server snapshots its browser assets at startup and sends `Cache-Control:
no-store`. A page cannot mix newly edited CSS with its earlier JavaScript from
the same running server. The build check reports a restart requirement when
application files on disk change. Save pending edits and finish active jobs
before restarting; the launcher never kills a job to obtain a newer interface.

This identity is not a comparison against GitHub, a dependency lock, or scientific
qualification. Those remain separate release checks. Existing run manifests and
original evidence are not rewritten by the launcher.

## Release practice

Keep the launcher, browser assets and Python service in the same reviewed
distribution. Test `workbench open` from a fresh installed wheel, check its About
record, and exercise a bundled offline example before approving an outgoing
archive. A development build identifier does not mean the changes have been
committed, synchronized to the publication checkout, pushed, or released.
