"""Parsing UniProt feature strings into catalytic and ligand-binding positions.

UniProt serves sequence features in a TSV column as a run-on string:

    ACT_SITE 195; /note="Charge relay system"; /evidence="ECO:0000255"; ACT_SITE 57; ...
    BINDING 57..59; /ligand="ATP"; /ligand_id="ChEBI:CHEBI:30616"

Everything downstream depends on reading these correctly, so the parser is
deliberately strict about two things:

* **A position it cannot read is dropped and counted, never guessed.** UniProt
  uses `?`, `<`, and `>` for uncertain and open-ended positions. A feature at an
  unknown position cannot be checked against a structure, and treating it as
  absent would understate the active site while treating it as present would
  overstate it. It is reported as unparsed instead.
* **Evidence codes are kept.** A catalytic residue asserted by similarity
  (ECO:0000250, ECO:0000256) is weaker evidence than one demonstrated
  experimentally (ECO:0000269), and the caller is entitled to know which it has.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Sequence

# Feature keys this module understands. Others in the column are ignored.
CATALYTIC_KEYS = ("ACT_SITE",)
LIGAND_KEYS = ("BINDING",)
OTHER_SITE_KEYS = ("SITE",)
KNOWN_KEYS = CATALYTIC_KEYS + LIGAND_KEYS + OTHER_SITE_KEYS

# Experimental evidence codes, as opposed to inferred-by-similarity ones.
EXPERIMENTAL_ECO = ("ECO:0000269", "ECO:0000305", "ECO:0007744")

_FEATURE_START = re.compile(r"\b(" + "|".join(KNOWN_KEYS) + r")\s+([^;]+)")
_QUALIFIER = re.compile(r'/(\w+)="([^"]*)"')
_POSITION = re.compile(r"^(\d+)(?:\.\.(\d+))?$")


@dataclass(frozen=True)
class Feature:
    """One annotated position or span."""

    key: str                      # ACT_SITE | BINDING | SITE
    start: int                    # 1-based, inclusive
    end: int                      # 1-based, inclusive; == start for a point feature
    note: str = ""
    ligand: str = ""
    ligand_id: str = ""
    evidence: str = ""

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(range(self.start, self.end + 1))

    @property
    def is_catalytic(self) -> bool:
        return self.key in CATALYTIC_KEYS

    @property
    def is_ligand(self) -> bool:
        return self.key in LIGAND_KEYS

    @property
    def experimentally_evidenced(self) -> bool:
        return any(code in self.evidence for code in EXPERIMENTAL_ECO)


@dataclass(frozen=True)
class FeatureSet:
    """Everything parsed for one protein, including what could not be read."""

    features: Sequence[Feature] = ()
    unparsed: Sequence[str] = field(default_factory=tuple)

    def of_kind(self, *keys: str) -> list[Feature]:
        return [f for f in self.features if f.key in keys]

    @property
    def catalytic(self) -> list[Feature]:
        return [f for f in self.features if f.is_catalytic]

    @property
    def ligand_binding(self) -> list[Feature]:
        return [f for f in self.features if f.is_ligand]

    def catalytic_positions(self) -> tuple[int, ...]:
        seen: list[int] = []
        for feature in self.catalytic:
            for position in feature.positions:
                if position not in seen:
                    seen.append(position)
        return tuple(sorted(seen))

    def ligand_positions(self) -> tuple[int, ...]:
        seen: list[int] = []
        for feature in self.ligand_binding:
            for position in feature.positions:
                if position not in seen:
                    seen.append(position)
        return tuple(sorted(seen))

    def ligands(self) -> tuple[str, ...]:
        names = [f.ligand for f in self.ligand_binding if f.ligand]
        return tuple(sorted(set(names)))

    @property
    def has_catalytic_annotation(self) -> bool:
        return bool(self.catalytic)


def _split_features(raw: str) -> list[tuple[str, str, int]]:
    """Split the run-on column into (key, body, offset) triples."""
    matches = list(_FEATURE_START.finditer(raw))
    out = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(raw)
        out.append((match.group(1), raw[match.start() : end], match.start()))
    return out


def parse_features(raw: str | None) -> FeatureSet:
    """Parse one UniProt feature column into positions and qualifiers."""
    if not raw or not raw.strip():
        return FeatureSet(features=(), unparsed=())

    features: list[Feature] = []
    unparsed: list[str] = []

    for key, body, _ in _split_features(raw):
        # The position token is whatever follows the key, up to the first ';'.
        after_key = body[len(key) :].lstrip()
        position_token = after_key.split(";", 1)[0].strip()

        match = _POSITION.match(position_token)
        if not match:
            # '?', '<12', '12..?' and similar: a position we cannot check.
            unparsed.append(f"{key} {position_token}")
            continue

        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else start
        if end < start:
            unparsed.append(f"{key} {position_token}")
            continue

        qualifiers = dict(_QUALIFIER.findall(body))
        features.append(
            Feature(
                key=key,
                start=start,
                end=end,
                note=qualifiers.get("note", ""),
                ligand=qualifiers.get("ligand", ""),
                ligand_id=qualifiers.get("ligand_id", ""),
                evidence=qualifiers.get("evidence", ""),
            )
        )

    return FeatureSet(features=tuple(features), unparsed=tuple(unparsed))


def parse_cofactors(raw: str | None) -> tuple[str, ...]:
    """Cofactor names from the UniProt `cc_cofactor` comment column.

    Shape: `COFACTOR: Name=Zn(2+); Xref=ChEBI:CHEBI:29105; Evidence=...;`
    """
    if not raw:
        return ()
    names = re.findall(r"Name=([^;]+);", raw)
    return tuple(sorted({name.strip() for name in names if name.strip()}))
