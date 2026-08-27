# Cross-platform testing and safe file sharing

The test bundle is designed for macOS and Linux with Python 3.10-3.12. It
contains synthetic examples and code; private campaign data and unpublished
sequences are not part of the bundle.

## Recipient steps

```bash
unzip yauvi-structural-workbench-prepublic-0.1.0-dev0.zip
cd yauvi-structural-workbench-prepublic-0.1.0-dev0
shasum -a 256 -c SHA256SUMS.txt       # macOS
# or: sha256sum -c SHA256SUMS.txt     # Linux

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[dev]"
python tools/run_structural_workbench_tests.py --core-only --json-out TEST_RESULT.json
```

Open `yauvi-structural-workbench/public-showcase/index.html` directly for the
offline narrative. The privacy-minimized test ZIP deliberately omits the larger
source-checkout controller and private-platform support surface. From the full
local source checkout, the loopback interface is:

```bash
structqc describe
```

## What the recipient should return

- `TEST_RESULT.json`
- operating system and version
- CPU architecture (`arm64`, `aarch64`, or `x86_64`)
- `python --version`
- installation errors, if any
- whether the static showcase opened
- whether the loopback interface started

Do not return private structures, sequences, patient data, credentials, or local
cache contents.

## File-sharing choices

For local testing, use a direct encrypted transfer, an access-controlled shared
folder, AirDrop, or removable media. Send the ZIP and its sidecar SHA-256 file as
separate items when practical, and ask the recipient to verify the checksum
before installation. Do not create a public repository, public link, archival
deposit, package-index release, or journal submission from this pre-public test
bundle.

The bundle is not an offline dependency mirror. The scientific examples run
offline after installation, but a new environment normally needs internet
access to obtain declared Python dependencies. Optional FreeSASA, Foldseek,
DIAMOND, MolProbity, Phenix, mkdssp, and MDAnalysis workflows require their own
installation and licenses.
