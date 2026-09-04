# Reproducing a panel on a second machine

Adoption protocol rule 7: a stratum passing only where its expectations were
recorded is not adopted. This is how someone else checks.

**You are not being asked to confirm the result.** A case that passes here and
fails on your machine is the finding, and it is worth more than agreement. The
membrane stratum passed 16/16 on the machine that recorded it and 9-10/16 on
every other machine tried, and that gap is the reason this step exists.

## What you need

- Python 3.11 or newer, with `numpy`, `scipy`, `gemmi` and `biopython` importable.
- For the sf-csa panel only: `foldseek` 10.941cd33 and `diamond` 2.1.11 on PATH.
  Other versions will produce different numbers rather than an error.
- No network. The panels run offline against checksummed files that ship with the
  repository, and `execution_policy.network_access` is `forbidden`.

## Running it

From this directory:

    sh reproduce_panel.sh ADOPTION_DRAFT_ABL.json
    sh reproduce_panel.sh ADOPTION_DRAFT_SFCSA.json

Each one builds its own virtual environment, installs the engine packages from
the source tree beside the distribution, verifies every artifact's SHA-256
against the panel's source lock, and executes. It prints where it wrote.

Expect the ABL panel to take under a minute. The sf-csa panel re-runs a full
twelve-query campaign per record and takes roughly ten.

## What to send back

    <output-directory>/EXECUTION_STATUS.json

Send the file, not a summary. It carries your platform, your Python version and
the measured numbers for every case, which is what identifies *which* case
diverged and by how much. A screenshot of the summary line loses exactly the
information that matters.

## What the output means

    14/14 cases passed, 1/1 controls passed      every case matched the recorded values
    coverage: 6/6 features witnessed             the panel exercised what it claims to
    scope_qualified: false                       always. This runner never qualifies a
                                                 scope; qualification needs this result
                                                 plus review, and is not something a
                                                 script decides.

A failing case prints the check that failed, what was expected and what was
observed. That is the useful part. Please send it even if — especially if — it
looks like your setup rather than the science: the recorded distances are
compared at a 0.001 A tolerance, and a genuine environment difference showing up
there is itself a result about how reproducible these numbers are.

## If it will not start

The script fails early and says why. The three common ones:

- *"could not install <package>"* — the engine sources are expected one level
  above the distribution directory. Check the repository was cloned whole.
- *"<cli> is not on PATH"* — an engine package installed but its console script
  did not. This has happened before and looked like a science failure.
- *"source lock failed"* — a file differs from its recorded digest. Do not work
  around it. That message is the check doing its job, and the mismatch is the
  thing to report.
