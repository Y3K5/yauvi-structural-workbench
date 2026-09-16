"""Fixtures for the activity-state tests."""
from __future__ import annotations

from pathlib import Path

import pytest

from actstate.core import ProteinRecord

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

# A canonical triad annotation, reused across tests.
TRIAD_ACT_SITE = (
    'ACT_SITE 5; /note="Charge relay system"; /evidence="ECO:0000269|PubMed:1"; '
    'ACT_SITE 9; /note="Charge relay system"; /evidence="ECO:0000255"; '
    'ACT_SITE 14; /note="Nucleophile"; /evidence="ECO:0000269|PubMed:1"'
)


def sequence_with(p5="H", p9="D", p14="S", length=20):
    seq = list("A" * length)
    seq[4], seq[8], seq[13] = p5, p9, p14
    return "".join(seq)


@pytest.fixture
def examples() -> Path:
    if not EXAMPLES.is_dir():
        pytest.skip("bundled examples are not present")
    return EXAMPLES


@pytest.fixture
def intact_record() -> ProteinRecord:
    return ProteinRecord(
        accession="P_TEST",
        sequence=sequence_with(),
        act_site_raw=TRIAD_ACT_SITE,
    )


@pytest.fixture
def write_pdb(tmp_path):
    """Build a small PDB from (seq_id, three_letter, (x, y, z)) tuples."""

    def build(name, residues, *, header="HEADER    TEST\n", heteroatoms=(), b=50.0):
        lines = [header]
        serial = 1
        for seq_id, res_name, (x, y, z) in residues:
            for atom in ("CA", "CB"):
                if atom == "CB" and res_name == "GLY":
                    continue
                offset = 0.0 if atom == "CA" else 0.5
                lines.append(
                    f"ATOM  {serial:>5} {atom:<4}{res_name:>4} A{seq_id:>4}    "
                    f"{x + offset:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{b:>6.2f}           C\n"
                )
                serial += 1
        for index, (res_name, (x, y, z), element) in enumerate(heteroatoms, start=1):
            lines.append(
                f"HETATM{serial:>5} {element:<4}{res_name:>4} A{900 + index:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 30.00          {element:>2}\n"
            )
            serial += 1
        lines.append("END\n")
        path = tmp_path / name
        path.write_text("".join(lines), encoding="utf-8")
        return path

    return build


@pytest.fixture
def clustered_triad():
    return [
        (5, "HIS", (0.0, 0.0, 0.0)),
        (9, "ASP", (4.0, 0.0, 0.0)),
        (14, "SER", (2.0, 3.5, 0.0)),
    ]


@pytest.fixture
def dispersed_triad():
    return [
        (5, "HIS", (0.0, 0.0, 0.0)),
        (9, "ASP", (30.0, 0.0, 0.0)),
        (14, "SER", (0.0, 30.0, 0.0)),
    ]
