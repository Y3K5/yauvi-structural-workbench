"""Classification of non-polymer components found in deposited structures.

Shipped as data, not hardcoded logic, so "what counts as a cryoprotectant" is
auditable and can be corrected without touching the code that uses it. Every
entry names the class it belongs to; anything absent is `unclassified` and is
never removed by default.

Ions are deliberately their own class. A magnesium may be a crystallisation
additive or the catalytic centre, and nothing in the three-letter code tells you
which — so ions are reported, not dropped, unless explicitly named.
"""
from __future__ import annotations

VERSION = "1.0"

SOLVENT = {"HOH", "DOD", "WAT", "D2O"}

CRYOPROTECTANT = {
    "GOL",  # glycerol
    "EDO",  # ethylene glycol
    "PEG", "PG4", "PGE", "P6G", "1PE", "2PE", "XPE",  # polyethylene glycols
    "MPD",  # 2-methyl-2,4-pentanediol
    "TRE", "SUC", "GLC",  # sugars used as cryoprotectant
    "DMS",  # dimethyl sulfoxide
}

BUFFER = {
    "TRS",  # tris
    "MES", "EPE", "HEPES", "BTB", "CIT", "TLA", "FLC",
    "ACT", "ACY",  # acetate
    "FMT",  # formate
    "IMD",  # imidazole
    "PO4", "SO4", "NO3", "CO3",  # commonly crystallisation salts
    "BME", "DTT",  # reducing agents
}

ION = {
    "NA", "K", "LI", "RB", "CS",
    "MG", "CA", "SR", "BA",
    "ZN", "FE", "FE2", "CU", "CU1", "MN", "CO", "NI", "CD", "HG",
    "CL", "BR", "IOD", "F",
}

CLASSES = {
    "solvent": SOLVENT,
    "cryoprotectant": CRYOPROTECTANT,
    "buffer": BUFFER,
    "ion": ION,
}

# Removed by default. Ions and unclassified components are reported, never
# dropped without being named: a structural metal and a crystallisation salt
# are indistinguishable by residue code alone.
DEFAULT_DROP = ("solvent", "cryoprotectant", "buffer")


def classify(resname: str) -> str:
    code = resname.strip().upper()
    for name, members in CLASSES.items():
        if code in members:
            return name
    return "unclassified"
