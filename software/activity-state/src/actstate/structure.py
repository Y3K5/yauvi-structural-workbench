"""A minimal structure reader: residues, coordinates, and heteroatoms.

Deliberately standard-library only. The house style for these packages is
`dependencies = []` — heavy tools are consumed through their outputs rather than
imported — and everything this module needs is a few fields at fixed columns.

It reads what the activity-state signals actually use and nothing more:

* per-residue identity and a representative side-chain coordinate
* heteroatom groups, separated into solvent, likely cryoprotectant/buffer, and
  everything else
* whether the coordinates are a prediction rather than an experiment

That last point is load-bearing. A predicted monomer is not evidence of a
functional conformation, and the classifier is required to refuse its strongest
label on prediction alone — so the reader must be able to say where the
coordinates came from.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Iterable, Sequence

# Ordered water and ion names that are not evidence of an occupied active site.
SOLVENT = frozenset({"HOH", "WAT", "DOD", "H2O"})

# Common crystallisation additives. Their presence is not cofactor occupancy,
# and counting them as such is the classic way to call an apo structure holo.
CRYO_AND_BUFFER = frozenset(
    {
        "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "MPD", "DMS",
        "SO4", "PO4", "CL", "NA", "K", "BR", "IOD", "ACT", "ACY", "FMT", "MES",
        "TRS", "EPE", "CIT", "TAR", "MLI", "NO3", "AZI", "CAC", "IMD", "BME",
    }
)

# Three-letter to one-letter, for checking an annotated residue's identity.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
}


class StructureError(RuntimeError):
    pass


@dataclass(frozen=True)
class Residue:
    """One amino acid, with the coordinate used for active-site geometry."""

    seq_id: int
    name: str
    chain: str
    x: float
    y: float
    z: float
    atom_used: str        # 'CB', or 'CA' for glycine and truncated side chains
    b_factor: float = 0.0

    @property
    def one_letter(self) -> str:
        return THREE_TO_ONE.get(self.name.upper(), "X")

    def distance_to(self, other: "Residue") -> float:
        return sqrt(
            (self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2
        )


@dataclass(frozen=True)
class Heteroatom:
    """A non-polymer group present in the coordinates."""

    name: str
    chain: str
    seq_id: int
    atom_count: int

    @property
    def is_solvent(self) -> bool:
        return self.name.upper() in SOLVENT

    @property
    def is_buffer_or_cryo(self) -> bool:
        return self.name.upper() in CRYO_AND_BUFFER

    @property
    def is_candidate_cofactor(self) -> bool:
        """Neither solvent nor a recognised additive — worth reporting."""
        return not self.is_solvent and not self.is_buffer_or_cryo


@dataclass(frozen=True)
class Structure:
    """Everything read from one coordinate file."""

    identifier: str
    residues: Sequence[Residue]
    heteroatoms: Sequence[Heteroatom]
    is_predicted: bool
    source_note: str = ""

    def by_seq_id(self, chain: str | None = None) -> dict[int, Residue]:
        return {
            r.seq_id: r for r in self.residues if chain is None or r.chain == chain
        }

    def chains(self) -> tuple[str, ...]:
        return tuple(sorted({r.chain for r in self.residues}))

    def candidate_cofactors(self) -> list[Heteroatom]:
        return [h for h in self.heteroatoms if h.is_candidate_cofactor]

    @property
    def mean_b_factor(self) -> float:
        if not self.residues:
            return 0.0
        return sum(r.b_factor for r in self.residues) / len(self.residues)


def _looks_predicted(text_head: str, path: Path) -> tuple[bool, str]:
    """Decide whether coordinates are a prediction, and say why."""
    upper = text_head.upper()
    if "ALPHAFOLD" in upper:
        return True, "header names AlphaFold"
    if "PREDICTED" in upper and "MODEL" in upper:
        return True, "header describes a predicted model"
    if path.name.upper().startswith("AF-"):
        return True, "filename follows the AlphaFold DB convention"
    if "BOLTZ" in upper or "ESMFOLD" in upper or "ROSETTAFOLD" in upper:
        return True, "header names a structure predictor"
    return False, ""


def _parse_pdb(text: str, *, identifier: str, path: Path) -> Structure:
    """Read ATOM/HETATM records at their fixed PDB columns."""
    # One representative atom per residue: CB where present, CA otherwise.
    # Glycine has no CB, and disordered side chains are often truncated to CA.
    chosen: dict[tuple[str, int], Residue] = {}
    het_atoms: dict[tuple[str, str, int], int] = {}

    for line in text.splitlines():
        record = line[:6]
        if record == "ENDMDL":
            break  # first model only; NMR ensembles are not averaged here
        if record not in ("ATOM  ", "HETATM"):
            continue
        try:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21:22].strip() or "A"
            seq_id = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except (ValueError, IndexError):
            continue
        try:
            b_factor = float(line[60:66])
        except (ValueError, IndexError):
            b_factor = 0.0

        if record == "HETATM" and res_name.upper() not in THREE_TO_ONE:
            key = (res_name, chain, seq_id)
            het_atoms[key] = het_atoms.get(key, 0) + 1
            continue

        if atom_name not in ("CB", "CA"):
            continue
        key = (chain, seq_id)
        existing = chosen.get(key)
        # CB wins; CA is only kept when no CB was seen.
        if existing is not None and existing.atom_used == "CB":
            continue
        chosen[key] = Residue(
            seq_id=seq_id, name=res_name, chain=chain,
            x=x, y=y, z=z, atom_used=atom_name, b_factor=b_factor,
        )

    predicted, note = _looks_predicted(text[:4000], path)
    return Structure(
        identifier=identifier,
        residues=tuple(sorted(chosen.values(), key=lambda r: (r.chain, r.seq_id))),
        heteroatoms=tuple(
            Heteroatom(name=name, chain=chain, seq_id=seq_id, atom_count=count)
            for (name, chain, seq_id), count in sorted(het_atoms.items())
        ),
        is_predicted=predicted,
        source_note=note,
    )


def _parse_mmcif(text: str, *, identifier: str, path: Path) -> Structure:
    """Read the atom_site loop of an mmCIF file.

    Only the columns needed here are resolved, by name rather than by position,
    because mmCIF does not fix column order.
    """
    lines = text.splitlines()
    headers: list[str] = []
    rows: list[str] = []
    in_loop = False
    collecting = False

    for line in lines:
        stripped = line.strip()
        if stripped == "loop_":
            in_loop, collecting, headers = True, False, []
            continue
        if in_loop and stripped.startswith("_atom_site."):
            headers.append(stripped.split(".", 1)[1])
            collecting = True
            continue
        if collecting:
            if not stripped or stripped.startswith(("#", "loop_", "_")):
                if rows:
                    break
                in_loop = collecting = False
                continue
            rows.append(stripped)

    if not headers or not rows:
        raise StructureError(f"no atom_site records found in {path}")

    index = {name: position for position, name in enumerate(headers)}
    required = ("group_PDB", "label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z")
    missing = [name for name in required if name not in index]
    if missing:
        raise StructureError(f"mmCIF atom_site is missing column(s): {', '.join(missing)}")

    res_name_col = index.get("label_comp_id", index.get("auth_comp_id"))
    chain_col = index.get("auth_asym_id", index.get("label_asym_id"))
    seq_col = index.get("auth_seq_id", index.get("label_seq_id"))
    b_col = index.get("B_iso_or_equiv")
    model_col = index.get("pdbx_PDB_model_num")

    chosen: dict[tuple[str, int], Residue] = {}
    het_atoms: dict[tuple[str, str, int], int] = {}
    first_model: str | None = None

    for row in rows:
        fields = row.split()
        if len(fields) < len(headers):
            continue
        if model_col is not None:
            if first_model is None:
                first_model = fields[model_col]
            elif fields[model_col] != first_model:
                break
        try:
            group = fields[index["group_PDB"]]
            atom_name = fields[index["label_atom_id"]].strip('"')
            res_name = fields[res_name_col] if res_name_col is not None else "UNK"
            chain = fields[chain_col] if chain_col is not None else "A"
            seq_raw = fields[seq_col] if seq_col is not None else "."
            if seq_raw in (".", "?"):
                continue
            seq_id = int(seq_raw)
            x = float(fields[index["Cartn_x"]])
            y = float(fields[index["Cartn_y"]])
            z = float(fields[index["Cartn_z"]])
        except (ValueError, IndexError, KeyError):
            continue
        b_factor = 0.0
        if b_col is not None:
            try:
                b_factor = float(fields[b_col])
            except (ValueError, IndexError):
                b_factor = 0.0

        if group == "HETATM" and res_name.upper() not in THREE_TO_ONE:
            key = (res_name, chain, seq_id)
            het_atoms[key] = het_atoms.get(key, 0) + 1
            continue
        if atom_name not in ("CB", "CA"):
            continue
        key = (chain, seq_id)
        existing = chosen.get(key)
        if existing is not None and existing.atom_used == "CB":
            continue
        chosen[key] = Residue(
            seq_id=seq_id, name=res_name, chain=chain,
            x=x, y=y, z=z, atom_used=atom_name, b_factor=b_factor,
        )

    predicted, note = _looks_predicted(text[:4000], path)
    return Structure(
        identifier=identifier,
        residues=tuple(sorted(chosen.values(), key=lambda r: (r.chain, r.seq_id))),
        heteroatoms=tuple(
            Heteroatom(name=name, chain=chain, seq_id=seq_id, atom_count=count)
            for (name, chain, seq_id), count in sorted(het_atoms.items())
        ),
        is_predicted=predicted,
        source_note=note,
    )


def read_structure(path: str | Path, *, identifier: str | None = None) -> Structure:
    """Read a PDB or mmCIF file. Format is chosen by content, then by suffix."""
    path = Path(path)
    if not path.is_file():
        raise StructureError(f"structure not found: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    name = identifier or path.stem

    head = text.lstrip()[:200]
    if head.startswith("data_") or path.suffix.lower() in (".cif", ".mmcif"):
        structure = _parse_mmcif(text, identifier=name, path=path)
    else:
        structure = _parse_pdb(text, identifier=name, path=path)

    if not structure.residues:
        raise StructureError(f"no amino-acid residues could be read from {path}")
    return structure


def pairwise_distances(residues: Iterable[Residue]) -> list[tuple[Residue, Residue, float]]:
    """Every pair and its separation, for cluster tests."""
    items = list(residues)
    return [
        (items[i], items[j], items[i].distance_to(items[j]))
        for i in range(len(items))
        for j in range(i + 1, len(items))
    ]
