"""The seam between measured evidence and the next workflow action.

Scientific modules measure.  Deterministic rules constrain.  A *decision provider*
chooses among the actions those rules left open.  This module defines that seam so
a future decision model can be plugged in without changing anything else, and
without ever being able to widen what it is allowed to choose.

The order is not negotiable:

    evidence summary -> hard rules -> permitted actions -> provider -> validation

A hard exclusion short-circuits the whole thing: when one has fired, no provider is
consulted at all, so a model cannot argue its way past a constraint it never sees.
Everything a provider returns is checked against the permitted set before it is
handed back, so a provider that misbehaves fails loudly rather than quietly
widening its own authority.

Core ships one provider: `DeterministicProvider`, which encodes no biology.  It
reports only whether the required evidence is present, absent, or in conflict.  A
learned provider (Laya-Y) is expected to be a second implementation of the same
protocol, registered under its own name.

A decision here is a *workflow action*, never a scientific conclusion.  "ADVANCE"
means the evidence permits the next analysis step; it does not mean the candidate
is protective, functional, correct, or real.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Protocol, Sequence, runtime_checkable
import math

#: Actions a provider may choose from.  Closed on purpose: a vocabulary that grows
#: silently makes two runs incomparable, and an unrecognised action cannot be
#: audited.  Extending it is a deliberate, versioned change.
DECISION_ACTIONS = (
    "ADVANCE",               # evidence permits the next approved step
    "DEPRIORITIZE",          # still plausible, weaker than alternatives
    "ADDITIONAL_ANALYSIS",   # a named analysis is missing and would resolve this
    "HUMAN_REVIEW",          # conflicting, unusual, or outside the provider's domain
    "INSUFFICIENT_EVIDENCE", # no justified conclusion either way, yet
)

#: Produced only by this module's rules, never by a provider.  Rejection is a
#: deterministic act; a learned model must not be able to discard a candidate.
HARD_EXCLUSION = "HARD_EXCLUSION"

#: Abstention is always available.  A provider that cannot justify a conclusion
#: must be able to say so, whatever else the rules permit.
ALWAYS_PERMITTED = ("HUMAN_REVIEW", "INSUFFICIENT_EVIDENCE")

#: How an evidence item stands.  "missing" and "conflicting" are first-class
#: states: silently treating either as a weak number is how a pipeline talks
#: itself into a conclusion it has not earned.
EVIDENCE_STRENGTHS = ("strong", "moderate", "weak", "conflicting", "missing")

#: Strengths that do not satisfy a requirement.
_UNSATISFIED = ("missing", "conflicting")

PROTOCOL_VERSION = "1.0"


class DecisionError(ValueError):
    """A state, provider, or returned decision violated the contract."""


@dataclass(frozen=True)
class ScientificState:
    """The compact summary a provider is allowed to see.

    Deliberately small.  A provider receives evidence *strengths*, not sequences,
    coordinates, or tool output, so that a decision cannot depend on anything the
    audit record does not also contain.
    """

    subject: str
    #: evidence name -> strength from EVIDENCE_STRENGTHS
    evidence: Mapping[str, str]
    #: evidence that must be present and unconflicted before ADVANCE is permitted
    required_evidence: Sequence[str] = ()
    #: hard exclusions that have already fired; any entry bypasses the provider
    exclusions: Sequence[str] = ()
    #: free-form provenance carried into the audit record (run id, workflow, …)
    context: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.subject or not str(self.subject).strip():
            raise DecisionError("a state must name its subject")
        if not self.evidence:
            raise DecisionError("a state must carry at least one evidence item")
        for name, strength in self.evidence.items():
            if not name or not str(name).strip():
                raise DecisionError("evidence names must be non-empty")
            if strength not in EVIDENCE_STRENGTHS:
                raise DecisionError(
                    f"unknown evidence strength {strength!r} for {name!r}; "
                    f"expected one of {', '.join(EVIDENCE_STRENGTHS)}")
        unknown = [name for name in self.required_evidence if name not in self.evidence]
        if unknown:
            # A requirement naming evidence the state never mentions is a bug in the
            # caller, not an absent measurement. Absent measurements are recorded as
            # "missing" so they appear in the audit record.
            raise DecisionError(
                "required evidence absent from the state entirely: " + ", ".join(sorted(unknown)))

    def unsatisfied_requirements(self) -> tuple[str, ...]:
        """Required evidence that is missing or in conflict."""
        return tuple(name for name in self.required_evidence
                     if self.evidence[name] in _UNSATISFIED)

    def conflicts(self) -> tuple[str, ...]:
        return tuple(sorted(n for n, s in self.evidence.items() if s == "conflicting"))

    def to_record(self) -> dict:
        return {
            "subject": self.subject,
            "evidence": dict(self.evidence),
            "required_evidence": list(self.required_evidence),
            "exclusions": list(self.exclusions),
            "context": dict(self.context),
        }


def permitted_actions(state: ScientificState) -> tuple[str, ...]:
    """The actions the deterministic rules leave open for this state.

    This is the constraint layer.  It is computed before any provider runs and is
    the only thing a provider may choose from.
    """
    if tuple(state.exclusions):
        # Nothing is open: the caller must return the hard exclusion instead.
        return ()
    allowed = [a for a in DECISION_ACTIONS if a in ALWAYS_PERMITTED]
    # A candidate may not advance while required evidence is missing or conflicting.
    # This is a rule, not a preference, which is why it lives here and not in a provider.
    if not state.unsatisfied_requirements():
        allowed.append("ADVANCE")
    allowed.append("DEPRIORITIZE")
    if state.unsatisfied_requirements():
        allowed.append("ADDITIONAL_ANALYSIS")
    return tuple(a for a in DECISION_ACTIONS if a in allowed)


@dataclass(frozen=True)
class ScientificDecision:
    """One typed decision, carrying enough to audit it later."""

    action: str
    rationale: str
    provider: str
    provider_version: str
    evidence_considered: Sequence[str] = ()
    missing_evidence: Sequence[str] = ()
    #: Optional action -> score. Labelled a preference distribution, not a
    #: probability: a number is only a probability once it has been calibrated
    #: against held-out outcomes, and nothing here has been.
    preferences: Mapping[str, float] | None = None

    def to_record(self) -> dict:
        record = {
            "protocol_version": PROTOCOL_VERSION,
            "action": self.action,
            "rationale": self.rationale,
            "provider": self.provider,
            "provider_version": self.provider_version,
            "evidence_considered": list(self.evidence_considered),
            "missing_evidence": list(self.missing_evidence),
            "interpretation": (
                "A workflow action selected from the permitted set. Not a scientific "
                "claim about the subject."),
        }
        if self.preferences is not None:
            record["preferences"] = dict(self.preferences)
            record["preferences_are_calibrated"] = False
        return record


@runtime_checkable
class DecisionProvider(Protocol):
    """What a decision provider must offer.

    A future Laya-Y adapter implements exactly this and nothing more: it receives a
    state and the permitted actions, and returns one of them with a reason.
    """

    provider_id: str
    version: str

    def decide(self, state: ScientificState, permitted: Sequence[str]) -> ScientificDecision:
        ...


class DeterministicProvider:
    """The reference provider: evidence bookkeeping, no biology, no model.

    It answers one question — does the evidence support taking the next step — using
    only presence, absence, and conflict.  It never returns a rejection: discarding a
    candidate is a deterministic act that belongs to the rules, not to a provider.
    """

    provider_id = "deterministic"
    version = "1.0"

    def decide(self, state: ScientificState, permitted: Sequence[str]) -> ScientificDecision:
        considered = tuple(sorted(state.evidence))
        unsatisfied = state.unsatisfied_requirements()
        conflicts = state.conflicts()

        def build(action: str, rationale: str) -> ScientificDecision:
            # Fall back to abstention when the rules did not leave the preferred
            # action open, rather than returning something impermissible.
            if action not in permitted:
                action, rationale = "INSUFFICIENT_EVIDENCE", (
                    f"{rationale} The preferred action was not permitted for this state.")
            return ScientificDecision(
                action=action, rationale=rationale,
                provider=self.provider_id, provider_version=self.version,
                evidence_considered=considered,
                missing_evidence=tuple(sorted(n for n, s in state.evidence.items()
                                              if s == "missing")))

        if conflicts:
            return build("HUMAN_REVIEW",
                         "Evidence disagrees and the disagreement is not resolvable by "
                         f"bookkeeping: {', '.join(conflicts)}.")
        if unsatisfied:
            return build("ADDITIONAL_ANALYSIS",
                         "Required evidence is not yet available: "
                         f"{', '.join(sorted(unsatisfied))}.")
        strengths = [state.evidence[name] for name in state.required_evidence]
        if strengths and all(s == "strong" for s in strengths):
            return build("ADVANCE",
                         "Every required evidence item is present and strong; no conflict "
                         "recorded. This permits the next step, and claims nothing further.")
        if any(s == "weak" for s in strengths):
            return build("DEPRIORITIZE",
                         "Required evidence is present but partly weak; other candidates "
                         "with stronger support should be preferred first.")
        return build("INSUFFICIENT_EVIDENCE",
                     "Evidence is present but not strong enough to justify a conclusion.")


_PROVIDERS: dict[str, DecisionProvider] = {}


def register_provider(provider: DecisionProvider, *, replace: bool = False) -> None:
    """Register a provider under its own id.

    Registration refuses to overwrite silently: two providers answering to one name
    would make an audit record ambiguous about which code produced a decision.
    """
    if not isinstance(provider, DecisionProvider):
        raise DecisionError("a provider must offer provider_id, version and decide()")
    name = provider.provider_id
    if not name or not str(name).strip():
        raise DecisionError("a provider must have a non-empty provider_id")
    if name in _PROVIDERS and not replace:
        raise DecisionError(f"a provider named {name!r} is already registered")
    _PROVIDERS[name] = provider


def get_provider(name: str) -> DecisionProvider:
    try:
        return _PROVIDERS[name]
    except KeyError:
        raise DecisionError(
            f"no decision provider named {name!r}; registered: "
            f"{', '.join(sorted(_PROVIDERS)) or 'none'}") from None


def registered_providers() -> tuple[str, ...]:
    return tuple(sorted(_PROVIDERS))


register_provider(DeterministicProvider())


def decide(state: ScientificState, provider: DecisionProvider | str = "deterministic") -> ScientificDecision:
    """Run the full sequence: rules, then provider, then validation.

    This is the only supported way to obtain a decision.  Calling a provider
    directly skips the constraint and validation steps that make its answer safe.
    """
    if isinstance(provider, str):
        provider = get_provider(provider)
    elif not isinstance(provider, DecisionProvider):
        raise DecisionError("a provider must offer provider_id, version and decide()")

    exclusions = tuple(state.exclusions)
    if exclusions:
        # The provider is never called. This is the point of the ordering.
        return ScientificDecision(
            action=HARD_EXCLUSION,
            rationale="Excluded by deterministic rule: " + "; ".join(exclusions),
            provider="core-rules", provider_version=PROTOCOL_VERSION,
            evidence_considered=tuple(sorted(state.evidence)))

    permitted = permitted_actions(state)
    decision = provider.decide(state, permitted)
    _validate(decision, state, permitted, provider)
    return decision


def _validate(decision: ScientificDecision, state: ScientificState,
              permitted: Sequence[str], provider: DecisionProvider) -> None:
    if not isinstance(decision, ScientificDecision):
        raise DecisionError("a provider must return a ScientificDecision")
    if decision.action == HARD_EXCLUSION:
        raise DecisionError("a provider may not issue a hard exclusion; that is a rule's decision")
    if decision.action not in DECISION_ACTIONS:
        raise DecisionError(f"unknown action {decision.action!r}")
    if decision.action not in permitted:
        raise DecisionError(
            f"provider {provider.provider_id!r} returned {decision.action!r}, which the rules "
            f"did not permit here; permitted: {', '.join(permitted)}")
    if not decision.rationale or not decision.rationale.strip():
        raise DecisionError("a decision must carry a rationale that can be read later")
    if decision.provider != provider.provider_id:
        raise DecisionError("a decision must be attributed to the provider that produced it")
    unknown = [name for name in decision.evidence_considered if name not in state.evidence]
    if unknown:
        # A provider citing evidence it was never given is either buggy or inventing
        # support for its answer; either way the audit record would be false.
        raise DecisionError("decision cites evidence absent from the state: " + ", ".join(sorted(unknown)))
    if decision.preferences is not None:
        bad = [a for a in decision.preferences if a not in permitted]
        if bad:
            raise DecisionError("preferences name actions that are not permitted: " + ", ".join(sorted(bad)))
        total = sum(decision.preferences.values())
        if any(v < 0 for v in decision.preferences.values()) or not math.isclose(total, 1.0, abs_tol=1e-6):
            raise DecisionError(f"preferences must be non-negative and sum to 1.0 (got {total})")


def describe() -> dict:
    """The declared contract for this seam, in the same shape modules use."""
    return {
        "module_id": "decision",
        "display_name": "Decision provider contract",
        "protocol_version": PROTOCOL_VERSION,
        "one_line": "Which workflow action do the rules and the evidence permit next?",
        "actions": list(DECISION_ACTIONS),
        "rule_only_actions": [HARD_EXCLUSION],
        "always_permitted": list(ALWAYS_PERMITTED),
        "evidence_strengths": list(EVIDENCE_STRENGTHS),
        "providers": list(registered_providers()),
        "inputs": [{"name": "state", "contract": "scientific_state",
                    "note": "evidence strengths, requirements, exclusions; no raw data"}],
        "outputs": [{"name": "decision", "contract": "scientific_decision", "format": "json"}],
        "limitations": [
            "A decision is a workflow action, not a scientific claim about the subject.",
            "A provider chooses only from actions the deterministic rules already permitted; "
            "it cannot widen that set.",
            "Hard exclusions bypass providers entirely, so a provider never sees a candidate "
            "it might argue should survive one.",
            "Rejection belongs to the rules. No provider may discard a candidate.",
            "The reference provider encodes no biology: it reports evidence completeness and "
            "conflict only.",
            "Preference scores are not calibrated probabilities and must not be read as "
            "reliability.",
        ],
    }
