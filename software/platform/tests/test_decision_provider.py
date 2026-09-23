"""A provider may choose among permitted actions and nothing else."""
from __future__ import annotations

import json

import pytest

from yauvi_platform.structural_workbench import decision as d


def state(**overrides):
    base = dict(
        subject="candidate-1",
        evidence={"structural_exposure": "strong", "conservation": "strong"},
        required_evidence=("structural_exposure", "conservation"),
    )
    base.update(overrides)
    return d.ScientificState(**base)


class Stub:
    """A provider that returns whatever the test tells it to."""

    provider_id = "stub"
    version = "0.1"

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.calls = 0

    def decide(self, st, permitted):
        self.calls += 1
        fields = dict(action="ADVANCE", rationale="because", provider=self.provider_id,
                      provider_version=self.version)
        fields.update(self.kwargs)
        return d.ScientificDecision(**fields)


# --- the state itself -------------------------------------------------------

def test_state_rejects_an_unknown_strength_and_an_empty_subject():
    with pytest.raises(d.DecisionError, match="unknown evidence strength"):
        d.ScientificState(subject="x", evidence={"a": "quite good"})
    with pytest.raises(d.DecisionError, match="name its subject"):
        d.ScientificState(subject="  ", evidence={"a": "strong"})
    with pytest.raises(d.DecisionError, match="at least one evidence item"):
        d.ScientificState(subject="x", evidence={})


def test_absent_evidence_is_recorded_as_missing_not_omitted():
    """Requiring evidence the state never mentions is a caller bug, not an absence."""
    with pytest.raises(d.DecisionError, match="absent from the state entirely"):
        d.ScientificState(subject="x", evidence={"a": "strong"}, required_evidence=("b",))
    st = state(evidence={"structural_exposure": "strong", "conservation": "missing"})
    assert st.unsatisfied_requirements() == ("conservation",)


# --- the constraint layer ---------------------------------------------------

def test_advance_is_not_permitted_while_a_requirement_is_unmet():
    ok = d.permitted_actions(state())
    assert "ADVANCE" in ok and "ADDITIONAL_ANALYSIS" not in ok
    blocked = d.permitted_actions(state(
        evidence={"structural_exposure": "strong", "conservation": "missing"}))
    assert "ADVANCE" not in blocked and "ADDITIONAL_ANALYSIS" in blocked


def test_abstention_is_always_available():
    for st in (state(), state(evidence={"structural_exposure": "conflicting",
                                        "conservation": "missing"})):
        assert set(d.ALWAYS_PERMITTED) <= set(d.permitted_actions(st))


def test_a_hard_exclusion_leaves_nothing_permitted_and_never_calls_the_provider():
    provider = Stub()
    excluded = state(exclusions=("host homology above threshold",))
    assert d.permitted_actions(excluded) == ()
    result = d.decide(excluded, provider)
    assert result.action == d.HARD_EXCLUSION
    assert result.provider == "core-rules"
    assert "host homology" in result.rationale
    assert provider.calls == 0, "a provider must not see an excluded candidate"


# --- validation of what a provider returns ----------------------------------

def test_a_provider_cannot_return_an_action_the_rules_withheld():
    unmet = state(evidence={"structural_exposure": "strong", "conservation": "missing"})
    with pytest.raises(d.DecisionError, match="did not permit"):
        d.decide(unmet, Stub(action="ADVANCE"))


def test_a_provider_cannot_issue_a_hard_exclusion_or_an_unknown_action():
    with pytest.raises(d.DecisionError, match="may not issue a hard exclusion"):
        d.decide(state(), Stub(action=d.HARD_EXCLUSION))
    with pytest.raises(d.DecisionError, match="unknown action"):
        d.decide(state(), Stub(action="DELETE_EVERYTHING"))


def test_a_decision_must_be_attributable_and_explained():
    with pytest.raises(d.DecisionError, match="rationale"):
        d.decide(state(), Stub(rationale="   "))
    with pytest.raises(d.DecisionError, match="attributed to the provider"):
        d.decide(state(), Stub(provider="someone-else"))


def test_a_provider_cannot_cite_evidence_it_was_never_given():
    with pytest.raises(d.DecisionError, match="cites evidence absent"):
        d.decide(state(), Stub(evidence_considered=("challenge_protection",)))


def test_preferences_must_cover_permitted_actions_and_sum_to_one():
    with pytest.raises(d.DecisionError, match="not permitted"):
        d.decide(state(), Stub(preferences={"ADDITIONAL_ANALYSIS": 1.0}))
    with pytest.raises(d.DecisionError, match="sum to 1.0"):
        d.decide(state(), Stub(preferences={"ADVANCE": 0.5, "HUMAN_REVIEW": 0.2}))
    good = d.decide(state(), Stub(preferences={"ADVANCE": 0.7, "HUMAN_REVIEW": 0.3}))
    assert good.to_record()["preferences_are_calibrated"] is False


def test_only_a_real_provider_is_accepted():
    with pytest.raises(d.DecisionError, match="provider_id"):
        d.decide(state(), object())


# --- the reference provider -------------------------------------------------

@pytest.mark.parametrize("evidence, expected", [
    ({"structural_exposure": "strong", "conservation": "strong"}, "ADVANCE"),
    ({"structural_exposure": "strong", "conservation": "missing"}, "ADDITIONAL_ANALYSIS"),
    ({"structural_exposure": "conflicting", "conservation": "strong"}, "HUMAN_REVIEW"),
    ({"structural_exposure": "strong", "conservation": "weak"}, "DEPRIORITIZE"),
    ({"structural_exposure": "moderate", "conservation": "moderate"}, "INSUFFICIENT_EVIDENCE"),
])
def test_reference_provider_reports_bookkeeping_not_biology(evidence, expected):
    result = d.decide(state(evidence=evidence))
    assert result.action == expected
    assert result.provider == "deterministic"
    assert result.rationale


def test_reference_provider_never_rejects_a_candidate():
    """Discarding a candidate is a rule's act, so no provider input should produce it."""
    for first in d.EVIDENCE_STRENGTHS:
        for second in d.EVIDENCE_STRENGTHS:
            result = d.decide(state(evidence={"structural_exposure": first,
                                              "conservation": second}))
            assert result.action in d.DECISION_ACTIONS
            assert result.action != d.HARD_EXCLUSION


def test_missing_evidence_is_named_in_the_decision():
    result = d.decide(state(evidence={"structural_exposure": "strong",
                                      "conservation": "missing"}))
    assert result.missing_evidence == ("conservation",)


# --- the registry and the audit record --------------------------------------

def test_registry_refuses_a_silent_overwrite_and_an_unknown_name():
    provider = Stub()
    d.register_provider(provider)
    try:
        with pytest.raises(d.DecisionError, match="already registered"):
            d.register_provider(Stub())
        d.register_provider(Stub(), replace=True)
        assert "stub" in d.registered_providers()
    finally:
        d._PROVIDERS.pop("stub", None)
    with pytest.raises(d.DecisionError, match="no decision provider named"):
        d.get_provider("laya-y")
    with pytest.raises(d.DecisionError, match="provider_id"):
        d.register_provider(object())


def test_the_audit_record_is_json_and_states_what_it_does_not_mean():
    record = d.decide(state()).to_record()
    text = json.dumps(record)
    assert json.loads(text)["action"] == "ADVANCE"
    assert record["protocol_version"] == d.PROTOCOL_VERSION
    assert "not a scientific claim" in record["interpretation"].lower()


def test_state_record_round_trips_for_provenance():
    record = state(context={"run_id": "RUN-1"}).to_record()
    assert json.loads(json.dumps(record))["context"]["run_id"] == "RUN-1"


def test_describe_declares_the_closed_vocabulary_and_its_limits():
    contract = describe = d.describe()
    assert contract["actions"] == list(d.DECISION_ACTIONS)
    assert contract["rule_only_actions"] == [d.HARD_EXCLUSION]
    assert "deterministic" in contract["providers"]
    assert any("not a scientific claim" in limit.lower() for limit in contract["limitations"])
    assert any("cannot widen" in limit.lower() for limit in contract["limitations"])
    json.dumps(describe)
