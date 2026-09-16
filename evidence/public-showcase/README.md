# YAUVI Structural Biology Platform — Mark 1 public evidence showcase

This directory is the public-safe, local microsite for **YAUVI Structural
Biology Platform — Mark 1**, powered by the YAUVI Structural Workbench. Open
`index.html` directly or use the loopback controller route `/public-showcase/`.

The HTML, CSS, and JavaScript are local. `data.js`, the sanitized evidence
copies, public-qualification summaries, reviewer pack, and checksum manifests
are generated from the five-case technical showcase, the reviewed SF-CSA
fixture case, recorded test baseline, release status, six analysis definitions,
and the checksum-locked qualification collection:

```bash
python tools/build_five_use_case_showcase.py --replace
python tools/verify_public_showcase.py
```

The first five demonstrations use synthetic scientific data. The synthetic
SF-CSA demonstration runs its canonical pipeline through deterministic Foldseek
and DIAMOND process-boundary test doubles that compute no alignments. A separate
qualification section reports the real local Foldseek and DIAMOND public
mini-case and the other five independent public cases. The two evidence classes
remain visibly separate. Building or viewing this directory does not publish,
deploy, upload, or submit anything.

The generated `share/` directory contains the canonical platform identity and
the start/share guide. These records keep the Mark 1 display name separate from
the unchanged standalone package and CLI names.
