# Method capabilities and readiness

Open the local Workbench, then choose **Methods · capability and readiness** in
the sidebar. No account, external model or hosted service is needed.

## What goes in and comes out

The view reads the existing method registry, named local source files and an
optional saved input case. It displays registration, source/configuration status,
input membership, the last recorded input check, expected outputs and the method's
claim ceiling. Selecting a case does not validate it or start a job. Use the
existing input check deliberately before execution.

Source presence, an executable path and a saved check are separate facts. This
view neither probes a runtime nor establishes dependency availability, tool network
behavior, execution success, scientific correctness or independent reproduction.
A changed source identity asks the user to refresh. Only named entry files are
fingerprinted; the card does not claim complete transitive source coverage.

## Distribution boundaries

The three cards are Structure QC, Collection research dossier and AntigenBench
synthetic check. The public distribution registers Structure QC. The other two
cards explicitly report **Unavailable in this distribution**: private collection
and AntigenBench adapters and their data are not shipped or substituted.

A broader local installation can report its already registered methods. The
AntigenBench card, where provided, describes an invented software example only;
it does not qualify real-protein evaluation or vaccine performance.

## Privacy and custody

Inspection uses local read-only API requests. It does not fetch public references,
import optional providers, probe engines, create examples, check inputs, queue jobs
or edit existing records. The view reports whether the host separately enables
public-reference retrieval. That setting is not changed by opening the view.

No new engine, provider, execution pathway, data store or cloud dependency is
introduced. Existing artifacts, validations, scientific qualification panels and
release gates retain their meanings.

## Check with invented inputs

`software/structural-workbench/tests/test_capability_distribution.py` checks absent
optional adapters, source drift, unchanged library bytes, forbidden execution and
network operations, guarded requests, incomplete input readiness and the served
native asset. `tools/test_workbench_wheel.py` independently rejects a wheel missing
the capability browser asset. These are software controls, not biological tests.
