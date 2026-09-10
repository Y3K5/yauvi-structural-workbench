"""structprep tests.

Two of these encode real mistakes that reached working code before this module
existed: a haem nearly stripped with the waters, and a biological assembly whose
duplicate chain ids made a tetramer look like a dimer.
"""
import json, os, pathlib, subprocess, sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from structprep.components import classify
from structprep.core import InputError, prepare, validate_policy, write_outputs

# Reference coordinates are not shipped: 3H8T is third-party RCSB data and this
# package makes no claim to redistribute it. Point STRUCTPREP_TEST_STRUCTURES at
# a directory holding 3H8T.pdb and 3H8T_bioassembly.pdb to run these; without it
# they skip. An absolute path was hardcoded here during development, which would
# have published a home directory and the name of an unrelated private study.
_ROOT = os.environ.get("STRUCTPREP_TEST_STRUCTURES")
CASE = pathlib.Path(_ROOT) if _ROOT else None
XRAY = CASE / "3H8T.pdb" if CASE else pathlib.Path("3H8T.pdb")
ASSEMBLY = CASE / "3H8T_bioassembly.pdb" if CASE else pathlib.Path("3H8T_bioassembly.pdb")

pytestmark = pytest.mark.skipif(
    not (CASE and XRAY.exists()),
    reason="set STRUCTPREP_TEST_STRUCTURES to a directory holding 3H8T.pdb",
)


def test_waters_go_and_the_cofactor_stays():
    """The whole point: strip solvent without stripping the functional ligand."""
    r = prepare(XRAY)
    kept = {a["resname"] for a in r["kept"] if a["record"] == "HETATM"}
    removed = {a["resname"] for a in r["removed"] if a["record"] == "HETATM"}
    assert "HEM" in kept, "the haem is the functional site; it must survive"
    assert "HOH" in removed
    assert r["counts"]["removed_by_class"].get("solvent", 0) > 400


def test_unclassified_components_are_never_dropped_by_default():
    """HEM is not in the table. Absence must mean 'keep', not 'unknown, discard'."""
    assert classify("HEM") == "unclassified"
    r = prepare(XRAY)
    assert not any(a["resname"] == "HEM" for a in r["removed"])


def test_contacting_additives_are_refused_not_removed():
    """Glycerol in a pocket may occupy a real site; removal must be deliberate."""
    r = prepare(XRAY)
    names = {x["component"] for x in r["refusals"]}
    assert {"GOL", "SO4"} & names, "additives touching the polymer should be refused"
    kept = {a["resname"] for a in r["kept"] if a["record"] == "HETATM"}
    assert "GOL" in kept


def test_explicit_drop_overrides_the_refusal():
    r = prepare(XRAY, {"drop": ["GOL", "SO4"]})
    kept = {a["resname"] for a in r["kept"] if a["record"] == "HETATM"}
    assert "GOL" not in kept and "SO4" not in kept
    assert "HEM" in kept


def test_solvent_is_exempt_from_contact_protection():
    """Nearly every crystal water touches protein; protecting them all would
    make the default policy a no-op."""
    r = prepare(XRAY)
    assert not any(x["class"] == "solvent" for x in r["refusals"])


def test_bridging_water_removal_is_warned():
    r = prepare(XRAY)
    assert any("bridges two chains" in w for w in r["warnings"])


def test_discarded_altloc_is_reported():
    r = prepare(XRAY)
    assert any("alternate location" in w for w in r["warnings"])


@pytest.mark.skipif(not ASSEMBLY.exists(), reason="assembly file absent")
def test_assembly_models_flatten_to_unique_chains(tmp_path):
    """The mistake this prevents: MODEL records reuse chain ids, so a tool
    reading chains alone sees one copy and reports crystal packing as if it
    were the biological interface."""
    r = prepare(ASSEMBLY)
    assert r["counts"]["models_in"] > 1
    write_outputs(r, tmp_path)
    rec = json.loads((tmp_path / "PREPARATION.json").read_text())
    assert len(rec["chain_remap"]) == 4
    chains = {l[21] for l in (tmp_path / "PREPARED.pdb").read_text().splitlines()
              if l.startswith("ATOM")}
    assert len(chains) == 4


def test_derivation_record_binds_parent_and_child(tmp_path):
    r = prepare(XRAY)
    write_outputs(r, tmp_path)
    rec = json.loads((tmp_path / "PREPARATION.json").read_text())
    assert rec["derived_from"]["sha256"] != rec["prepared_sha256"]
    assert len(rec["derived_from"]["sha256"]) == 64
    # the record must say the output is not docking-ready
    assert any("protonation" in l for l in rec["limitations"])


def test_unknown_class_is_rejected():
    with pytest.raises(InputError):
        validate_policy({"drop_classes": ["solvent", "nonsense"]})


def test_run_exits_nonzero_when_something_was_refused(tmp_path):
    """A caller who asked for a removal that did not happen should learn it from
    the exit code, not by reading a file."""
    env = {"PYTHONPATH": str(pathlib.Path(__file__).resolve().parents[1] / "src")}
    p = subprocess.run([sys.executable, "-m", "structprep.cli", "run",
                        "--structure", str(XRAY), "--out", str(tmp_path)],
                       capture_output=True, text=True, env={**env})
    assert p.returncode == 1
    assert json.loads(p.stdout)["refusals"] > 0
