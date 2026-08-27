"""Does the available evidence support this protein being in a working state?

Five signals are computed and **reported separately**. They are never summed,
averaged, or collapsed into a single score — the same rule `sf_csa` follows for
structural versus sequence similarity, and for the same reason: the signals
answer different questions, fail in different ways, and a reader who is told only
the total cannot tell which one carried it.

The label vocabulary is closed. Six values, and no free text is ever promoted
into one:

    active_state_supported     every signal that could be evaluated is
                               consistent with a working state, on experimental
                               coordinates, and none is unavailable. Note the
                               exact claim: a signal may be `unevaluated`
                               (the question does not arise for this entry —
                               no declared cofactor, or a single annotated
                               position with no geometry to check) and the
                               label still stands. It is not a guarantee that
                               all five signals fired.
    probable_active            residues intact, but a signal that would be needed
                               to say more is unavailable or predicted
    apo_but_competent          the site is intact and unoccupied; a cofactor the
                               entry declares is not present in the coordinates
    inactive_conformation      residues present but not mutually positioned as a
                               site
    active_site_disrupted      an annotated catalytic position does not hold a
                               residue that can perform chemistry
    indeterminate              not enough is annotated to make any claim

Three rules constrain what may be concluded:

1. **No annotated active site means `indeterminate`, never `inactive`.** Absence
   of annotation is not evidence of absence of function; most proteins in a
   proteome have no ACT_SITE line and are not thereby pseudoenzymes.
2. **A predicted model alone can never yield `active_state_supported`.** A
   predictor reproduces the fold it was trained to reproduce; it is not an
   observation of a functional state.
3. **An unavailable signal is recorded as unavailable.** It never becomes a
   neutral or favourable value, and it never silently drops out of the summary.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

from .features import FeatureSet, parse_cofactors, parse_features
from .structure import Structure, pairwise_distances

# Residues that can participate in catalysis — nucleophiles, acid/base pairs,
# metal ligands, and the two that act through backbone geometry. A position
# annotated ACT_SITE holding anything outside this set is the signature of a
# degraded site, which is what pseudoenzymes look like.
CATALYTICALLY_COMPETENT = frozenset("DEHCKRSTYNQGW")

# Upper bound on the separation of two residues that belong to one active site.
# Generous on purpose: catalytic residues are typically within ~10 A, and the
# question here is whether they are clustered at all, not how tightly.
SITE_CLUSTER_MAX_ANGSTROM = 16.0

LABELS = (
    "active_state_supported",
    "probable_active",
    "apo_but_competent",
    "inactive_conformation",
    "active_site_disrupted",
    "indeterminate",
)

# What a signal can say about itself.
SIGNAL_STATES = ("supported", "contradicted", "unevaluated", "unavailable")


@dataclass(frozen=True)
class Signal:
    """One line of evidence, with its own state and its own reason."""

    name: str
    state: str
    detail: str
    values: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.state not in SIGNAL_STATES:
            raise ValueError(f"unknown signal state: {self.state}")
        if not self.detail:
            raise ValueError(f"signal {self.name} must explain itself")

    @property
    def informative(self) -> bool:
        return self.state in ("supported", "contradicted")


@dataclass(frozen=True)
class ProteinRecord:
    """One protein's annotation, as read from a UniProt-style export."""

    accession: str
    sequence: str = ""
    act_site_raw: str = ""
    binding_raw: str = ""
    site_raw: str = ""
    cofactor_raw: str = ""
    ec_number: str = ""
    interpro: str = ""
    pfam: str = ""
    protein_name: str = ""

    def features(self) -> FeatureSet:
        joined = "; ".join(p for p in (self.act_site_raw, self.binding_raw, self.site_raw) if p)
        return parse_features(joined)

    def declared_cofactors(self) -> tuple[str, ...]:
        return parse_cofactors(self.cofactor_raw)


@dataclass(frozen=True)
class ActivityAssessment:
    """The result for one protein: a label, and the signals that produced it."""

    accession: str
    label: str
    signals: Sequence[Signal]
    catalytic_positions: Sequence[int] = ()
    declared_cofactors: Sequence[str] = ()
    unparsed_features: Sequence[str] = ()
    rationale: str = ""

    def __post_init__(self) -> None:
        if self.label not in LABELS:
            raise ValueError(f"label outside the closed vocabulary: {self.label}")

    def signal(self, name: str) -> Signal | None:
        return next((s for s in self.signals if s.name == name), None)


# -- the five signals -----------------------------------------------------


def completeness_signal(record: ProteinRecord, features: FeatureSet) -> Signal:
    """Are the annotated catalytic positions present, and can they do chemistry?"""
    positions = features.catalytic_positions()
    if not positions:
        return Signal(
            "completeness",
            "unevaluated",
            "no ACT_SITE annotation for this entry, so site integrity cannot be assessed",
        )
    if not record.sequence:
        return Signal(
            "completeness",
            "unavailable",
            f"{len(positions)} catalytic position(s) annotated but no sequence was supplied",
            {"catalytic_positions": list(positions)},
        )

    length = len(record.sequence)
    out_of_range = [p for p in positions if p < 1 or p > length]
    if out_of_range:
        return Signal(
            "completeness",
            "unavailable",
            (
                f"annotated position(s) {out_of_range} fall outside the supplied sequence "
                f"of length {length}; the sequence and the annotation are not the same entry"
            ),
            {"out_of_range": out_of_range, "sequence_length": length},
        )

    observed = {p: record.sequence[p - 1].upper() for p in positions}
    degraded = {p: aa for p, aa in observed.items() if aa not in CATALYTICALLY_COMPETENT}
    if degraded:
        return Signal(
            "completeness",
            "contradicted",
            (
                "annotated catalytic position(s) hold residues that cannot perform "
                "chemistry: "
                + ", ".join(f"{p}{aa}" for p, aa in sorted(degraded.items()))
            ),
            {"degraded": {str(p): aa for p, aa in degraded.items()}, "observed": {str(p): aa for p, aa in observed.items()}},
        )

    experimental = sum(1 for f in features.catalytic if f.experimentally_evidenced)
    return Signal(
        "completeness",
        "supported",
        (
            f"all {len(positions)} annotated catalytic position(s) hold competent residues "
            f"({', '.join(f'{p}{aa}' for p, aa in sorted(observed.items()))}); "
            f"{experimental} of {len(features.catalytic)} annotation(s) carry experimental evidence"
        ),
        {
            "observed": {str(p): aa for p, aa in observed.items()},
            "experimentally_evidenced": experimental,
        },
    )


def geometry_signal(
    features: FeatureSet,
    structure: Structure | None,
    *,
    chain: str | None = None,
    max_separation: float = SITE_CLUSTER_MAX_ANGSTROM,
) -> Signal:
    """Are the catalytic residues mutually positioned as one site?"""
    positions = features.catalytic_positions()
    if not positions:
        return Signal("geometry", "unevaluated", "no catalytic positions to place")
    if structure is None:
        return Signal(
            "geometry",
            "unavailable",
            "no structure supplied, so active-site geometry was not evaluated",
        )
    if len(positions) < 2:
        return Signal(
            "geometry",
            "unevaluated",
            "a single catalytic position has no geometry to check against itself",
            {"catalytic_positions": list(positions)},
        )

    # A structure with several chains carrying the same annotated positions is
    # ambiguous, and `by_seq_id` keys on seq_id alone: without a chain, one
    # chain's residues silently overwrite another's and the geometry below would
    # describe whichever chain was parsed last. Say so instead of guessing.
    if chain is None:
        carrying = sorted(
            {r.chain for r in structure.residues if r.seq_id in set(positions)}
        )
        if len(carrying) > 1:
            return Signal(
                "geometry",
                "unavailable",
                (
                    f"catalytic positions appear in {len(carrying)} chains "
                    f"({', '.join(carrying)}) and no chain was selected; geometry "
                    "would be measured on an arbitrary one, so it was not evaluated"
                ),
                {"chains_carrying_catalytic_positions": carrying},
            )

    residues = structure.by_seq_id(chain)
    found = [residues[p] for p in positions if p in residues]
    missing = [p for p in positions if p not in residues]
    if len(found) < 2:
        return Signal(
            "geometry",
            "unavailable",
            (
                f"only {len(found)} of {len(positions)} catalytic position(s) are resolved in "
                f"the coordinates (missing {missing}); geometry cannot be evaluated"
            ),
            {"missing_from_structure": missing},
        )

    pairs = pairwise_distances(found)
    widest = max(pairs, key=lambda item: item[2])
    max_distance = round(widest[2], 2)
    detail_positions = ", ".join(str(r.seq_id) for r in found)

    if max_distance > max_separation:
        return Signal(
            "geometry",
            "contradicted",
            (
                f"catalytic residues ({detail_positions}) are not clustered: widest separation "
                f"{max_distance} A between {widest[0].seq_id} and {widest[1].seq_id}, "
                f"above the {max_separation} A bound for one site"
            ),
            {"max_separation_angstrom": max_distance, "resolved": len(found), "missing_from_structure": missing},
        )
    return Signal(
        "geometry",
        "supported",
        (
            f"catalytic residues ({detail_positions}) form a cluster: widest separation "
            f"{max_distance} A"
            + (f"; {len(missing)} position(s) unresolved in the coordinates" if missing else "")
        ),
        {"max_separation_angstrom": max_distance, "resolved": len(found), "missing_from_structure": missing},
    )


def occupancy_signal(record: ProteinRecord, structure: Structure | None) -> Signal:
    """Apo or holo: does the entry declare a cofactor, and is one present?"""
    declared = record.declared_cofactors()
    if structure is None:
        if declared:
            return Signal(
                "occupancy",
                "unavailable",
                f"entry declares cofactor(s) {list(declared)} but no structure was supplied",
                {"declared_cofactors": list(declared)},
            )
        return Signal("occupancy", "unevaluated", "no declared cofactor and no structure supplied")

    present = structure.candidate_cofactors()
    names = sorted({h.name for h in present})

    if not declared and not present:
        return Signal(
            "occupancy",
            "unevaluated",
            "no cofactor is declared and none is present; occupancy is not a question here",
        )
    if declared and not present:
        return Signal(
            "occupancy",
            "contradicted",
            (
                f"entry declares cofactor(s) {list(declared)} but the coordinates contain no "
                f"non-solvent heteroatom group; this is an apo structure"
            ),
            {"declared_cofactors": list(declared), "present": []},
        )
    if present and not declared:
        return Signal(
            "occupancy",
            "supported",
            (
                f"coordinates contain heteroatom group(s) {names} that are neither solvent nor a "
                f"recognised buffer or cryoprotectant, though the entry declares no cofactor"
            ),
            {"declared_cofactors": [], "present": names},
        )
    return Signal(
        "occupancy",
        "supported",
        f"declared cofactor(s) {list(declared)} and heteroatom group(s) {names} are both present",
        {"declared_cofactors": list(declared), "present": names},
    )


def conformation_signal(comparison: Mapping[str, object] | None) -> Signal:
    """How do these coordinates compare to references of known state?

    The comparison is produced by a structural aligner (Foldseek), which is an
    external runtime. When it has not been run the signal is `unavailable` — the
    fail-closed behaviour `shared/runtime-registry.yaml` mandates. It never
    becomes a neutral value that quietly stops mattering.
    """
    if not comparison:
        return Signal(
            "conformation",
            "unavailable",
            (
                "no reference-state comparison was supplied; a structural aligner and a "
                "curated set of references of known state are required to evaluate this"
            ),
        )
    reference = str(comparison.get("reference", "")).strip()
    state = str(comparison.get("state", "")).strip()
    if not reference or state not in ("active", "inactive"):
        return Signal(
            "conformation",
            "unavailable",
            (
                "reference comparison is present but does not name a reference and a state of "
                "'active' or 'inactive'; it cannot be interpreted"
            ),
            dict(comparison),
        )
    values = {"reference": reference, "reference_state": state}
    if "score" in comparison:
        values["score"] = comparison["score"]
    if state == "active":
        return Signal(
            "conformation",
            "supported",
            f"closest reference of known state is {reference}, which is an active-state structure",
            values,
        )
    return Signal(
        "conformation",
        "contradicted",
        f"closest reference of known state is {reference}, which is an inactive-state structure",
        values,
    )


def assembly_signal(fold_state: Mapping[str, object] | None) -> Signal:
    """Does the evidence describe the isolated fold or the working assembly?

    Composes with the platform's existing `fold_state` module rather than
    duplicating it: that module says *which* fold the evidence describes, this
    one asks whether that fold is competent.
    """
    if not fold_state:
        return Signal(
            "assembly",
            "unavailable",
            (
                "no fold_state record supplied; whether these coordinates describe the isolated "
                "monomer or the working assembly is undetermined"
            ),
        )
    state = str(fold_state.get("state", "")).strip()
    if state == "active_assembly":
        return Signal(
            "assembly",
            "supported",
            "fold_state reports the evidence describes the protein in its working assembly",
            dict(fold_state),
        )
    if state == "isolated_fold":
        return Signal(
            "assembly",
            "contradicted",
            (
                "fold_state reports the evidence describes the isolated monomer fold, not the "
                "assembly the protein works in"
            ),
            dict(fold_state),
        )
    return Signal(
        "assembly",
        "unavailable",
        f"fold_state record does not name a recognised state (got {state!r})",
        dict(fold_state),
    )


# -- combination ----------------------------------------------------------


def assign_label(signals: Sequence[Signal], *, structure: Structure | None) -> tuple[str, str]:
    """Choose one label from the closed vocabulary, and say why.

    Ordered by what each signal can rule out. A contradicted completeness signal
    is decisive on its own: if a catalytic position cannot do chemistry, nothing
    downstream can restore it.
    """
    by_name = {s.name: s for s in signals}
    completeness = by_name["completeness"]
    geometry = by_name["geometry"]
    occupancy = by_name["occupancy"]
    conformation = by_name["conformation"]
    assembly = by_name["assembly"]

    if completeness.state == "contradicted":
        return "active_site_disrupted", completeness.detail

    if completeness.state == "unevaluated":
        # Rule 1: nothing annotated means nothing to conclude.
        return (
            "indeterminate",
            "no catalytic site is annotated for this entry, so no activity-state claim is made",
        )

    if completeness.state == "unavailable":
        return "indeterminate", completeness.detail

    # From here, completeness is supported.
    if geometry.state == "contradicted":
        return "inactive_conformation", geometry.detail
    if conformation.state == "contradicted":
        return "inactive_conformation", conformation.detail

    if occupancy.state == "contradicted":
        return "apo_but_competent", occupancy.detail

    # Everything from here caps the claim at `probable_active` rather than
    # inverting it, and every reason for the cap is named in the rationale. They
    # are collected rather than returned one at a time: a rationale that names
    # only the first reason it hit reads as though the others were checked and
    # passed.
    reasons: list[str] = []

    if structure is not None and structure.is_predicted:
        # Rule 2: a prediction is not an observation of a functional state.
        # Collected rather than returned immediately: a predicted model whose
        # fold_state ALSO says isolated monomer has two reasons for the cap, and
        # returning here named only the first.
        reasons.append(
            "the coordinates are a predicted model"
            + (f" ({structure.source_note})" if structure.source_note else "")
            + ", which cannot on its own establish a functional state"
        )

    if assembly.state == "contradicted":
        # The evidence describes the isolated monomer, not the assembly the
        # protein works in. That is a limit on what was observed rather than a
        # finding about the protein, so it caps the claim instead of inverting
        # it: `inactive_conformation` would assert the monomer fold is a
        # non-functional pose, which fold_state does not say.
        reasons.append(assembly.detail)

    unavailable = sorted(s.name for s in signals if s.state == "unavailable")
    if unavailable:
        reasons.append(", ".join(unavailable) + " could not be evaluated")

    if reasons:
        return (
            "probable_active",
            (
                "the site is intact, but "
                + "; and ".join(reasons)
                + ", so a working state is not established"
            ),
        )

    return (
        "active_state_supported",
        "every evaluated signal is consistent with a working state on experimental coordinates",
    )


def assess(
    record: ProteinRecord,
    *,
    structure: Structure | None = None,
    chain: str | None = None,
    reference_comparison: Mapping[str, object] | None = None,
    fold_state: Mapping[str, object] | None = None,
    max_separation: float = SITE_CLUSTER_MAX_ANGSTROM,
) -> ActivityAssessment:
    """Assess one protein. Pure: no IO, no network, no paths."""
    features = record.features()
    signals = (
        completeness_signal(record, features),
        geometry_signal(features, structure, chain=chain, max_separation=max_separation),
        occupancy_signal(record, structure),
        conformation_signal(reference_comparison),
        assembly_signal(fold_state),
    )
    label, rationale = assign_label(signals, structure=structure)
    return ActivityAssessment(
        accession=record.accession,
        label=label,
        signals=signals,
        catalytic_positions=features.catalytic_positions(),
        declared_cofactors=record.declared_cofactors(),
        unparsed_features=features.unparsed,
        rationale=rationale,
    )
