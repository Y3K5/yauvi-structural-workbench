# The decision provider contract

The workbench's modules answer measurement questions: how trustworthy are these
coordinates, which residues sit in the membrane, which residues form an interface.
Something still has to answer the question that follows — **what happens to this
candidate next?**

This page describes the seam where that answer is produced, why it is shaped the way
it is, and how to supply your own.

## The sequence

```
evidence summary  ->  hard rules  ->  permitted actions  ->  provider  ->  validation
```

Each step exists to constrain the next one.

1. **The evidence summary** is small on purpose. A provider sees evidence *strengths*
   — `strong`, `moderate`, `weak`, `conflicting`, `missing` — not sequences,
   coordinates or tool output. A decision therefore cannot depend on anything the
   audit record does not also contain.
2. **Hard rules** run first and compute the permitted actions. They are deterministic
   and are not supplied by the provider.
3. **The provider chooses** one permitted action and states why.
4. **Validation** checks the answer against the permitted set before it is returned.

## What a provider may and may not do

| | |
|---|---|
| May | Choose among the actions the rules left open |
| May | Abstain — `HUMAN_REVIEW` and `INSUFFICIENT_EVIDENCE` are always available |
| May | Report preference scores across permitted actions |
| **May not** | Return an action the rules withheld |
| **May not** | Reject a candidate — rejection is a rule's act, never a provider's |
| **May not** | See a candidate that a hard exclusion already removed |
| **May not** | Cite evidence it was not given |

The last four are enforced, not advised. Each one has a test.

### Hard exclusions bypass providers entirely

When a deterministic exclusion has fired, the permitted set is empty and the provider
is **never called**. A model cannot argue its way past a constraint it never sees.
This is the single most important property of the design, and the reason the ordering
above is fixed.

## The vocabulary

Actions a provider may choose from:

| Action | Meaning |
|---|---|
| `ADVANCE` | The evidence permits the next approved step |
| `DEPRIORITIZE` | Still plausible, weaker than alternatives |
| `ADDITIONAL_ANALYSIS` | A named analysis is missing and would resolve this |
| `HUMAN_REVIEW` | Conflicting, unusual, or outside the provider's domain |
| `INSUFFICIENT_EVIDENCE` | No justified conclusion either way, yet |

`HARD_EXCLUSION` also exists, and is produced only by the rules.

The vocabulary is closed. A set that grows silently makes two runs incomparable, and
an unrecognised action cannot be audited.

**`ADVANCE` is a workflow action, not a scientific claim.** It means the evidence
permits the next analysis step. It does not mean the candidate is functional,
protective, correct or real.

## The reference provider

The package ships one provider, `deterministic`. It encodes **no biology**. It reports
only whether required evidence is present, absent, or in conflict:

| Evidence state | Decision |
|---|---|
| Anything conflicting | `HUMAN_REVIEW` |
| A requirement missing | `ADDITIONAL_ANALYSIS`, naming what is missing |
| All requirements strong | `ADVANCE` |
| A requirement weak | `DEPRIORITIZE` |
| Otherwise | `INSUFFICIENT_EVIDENCE` |

It is deliberately dull. Its job is to make the seam usable and testable without
asserting anything about biology.

## Using it

```python
from yauvi_platform.structural_workbench import ScientificState, decide

state = ScientificState(
    subject="candidate-1",
    evidence={"structural_exposure": "strong", "conservation": "missing"},
    required_evidence=("structural_exposure", "conservation"),
    context={"run_id": "RUN-2026-000143"},
)

result = decide(state)             # -> ADDITIONAL_ANALYSIS
result.to_record()                 # JSON-ready audit record
```

Evidence that was not measured is recorded as `missing`, never omitted. An omission
disappears from the audit record; a `missing` entry stays visible and is named in the
decision.

## Supplying your own provider

A provider needs three things: an id, a version, and a `decide` method that takes the
state and the permitted actions and returns a decision.

```python
class MyProvider:
    provider_id = "my-provider"
    version = "0.1"

    def decide(self, state, permitted):
        ...  # must return a ScientificDecision whose action is in `permitted`

register_provider(MyProvider())
decide(state, "my-provider")
```

Always go through `decide()`. Calling a provider directly skips the constraint and
validation steps that make its answer safe.

Registration refuses to overwrite an existing name unless you pass `replace=True`,
because two providers answering to one name would make an audit record ambiguous about
which code produced a decision.

## Preference scores are not probabilities

A provider may return scores across the permitted actions. They must be non-negative
and sum to 1. The audit record labels them `preferences_are_calibrated: false`, and
that stays false until scores have been checked against held-out outcomes. A number is
only a probability once it has been calibrated; until then it is a model preference and
is labelled as one.

## Why this exists now

This contract is an extension point, not a feature. It lets a future learned decision
layer be added as a second provider without changing the packages around it — and,
just as importantly, it fixes *now*, while the surface is small and public, exactly
how much authority such a layer is ever allowed to have.
