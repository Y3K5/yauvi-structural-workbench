"""The accessibility palette must stay separable for dichromatic viewers.

The figure this palette colours is the one a reader decides from: which residues
are antibody-accessible, and which are buried. If two classes collapse for a
reader with colour vision deficiency, the figure misleads exactly the person
using it most carefully.

These tests re-derive the choice rather than pin the hex values, so a future
change is measured against the same criterion instead of eyeballed. An earlier
by-eye assignment of the same palette scored *worse* than what it replaced.
"""
from __future__ import annotations

import itertools
import math
import pathlib

import pytest

from memorient.viz import ACC_COLORS, ACC_COLORS_LEGACY, apply_current_palette


# --- Vienot, Brettel & Mollon (1999) dichromat simulation ------------------

_M = [[17.8824, 43.5161, 4.11935],
      [3.45565, 27.1554, 3.86714],
      [0.0299566, 0.184309, 1.46709]]
_MINV = [[0.080944, -0.130504, 0.116721],
         [-0.0102485, 0.0540194, -0.113615],
         [-0.000365294, -0.00412163, 0.693513]]
_SIM = {
    "normal": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "deuteranope": [[1, 0, 0], [0.494207, 0, 1.24827], [0, 0, 1]],
    "protanope": [[0, 2.02344, -2.52581], [0, 1, 0], [0, 0, 1]],
    "tritanope": [[1, 0, 0], [0, 1, 0], [-0.395913, 0.801109, 0]],
}


def _srgb_to_linear(value: float) -> float:
    return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4


def _linear_to_srgb(value: float) -> float:
    value = max(0.0, min(1.0, value))
    return 12.92 * value if value <= 0.0031308 else 1.055 * value ** (1 / 2.4) - 0.055


def _mul(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def _simulate(hex_color: str, kind: str):
    rgb = [int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [_srgb_to_linear(c) for c in rgb]
    out = _mul(_MINV, _mul(_SIM[kind], _mul(_M, linear)))
    return tuple(_linear_to_srgb(c) for c in out)


def _lab(rgb):
    linear = [_srgb_to_linear(c) for c in rgb]
    x = 0.4124 * linear[0] + 0.3576 * linear[1] + 0.1805 * linear[2]
    y = 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
    z = 0.0193 * linear[0] + 0.1192 * linear[1] + 0.9505 * linear[2]

    def f(t):
        return t ** (1 / 3) if t > (6 / 29) ** 3 else t / (3 * (6 / 29) ** 2) + 4 / 29

    fx, fy, fz = f(x / 0.95047), f(y / 1.0), f(z / 1.08883)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))


def worst_case_separation(colors) -> float:
    """Smallest perceptual distance between any two colours, over all vision types."""
    return min(
        math.dist(_lab(_simulate(a, kind)), _lab(_simulate(b, kind)))
        for kind in _SIM
        for a, b in itertools.combinations(colors, 2)
    )


# --- the properties that must hold -----------------------------------------

#: Below this, two classes are close enough that a reader could confuse them.
#: CIE76 dE of 10 is the conventional "clearly different" threshold.
MINIMUM_SEPARATION = 10.0


def test_every_class_is_separable_for_every_vision_type():
    separation = worst_case_separation(list(ACC_COLORS.values()))
    assert separation >= MINIMUM_SEPARATION, (
        f"worst-case separation {separation:.1f} is below {MINIMUM_SEPARATION}; "
        "two accessibility classes are confusable"
    )


def test_current_palette_beats_the_one_it_replaced():
    """The change has to be an improvement, not merely a different set of hues."""
    assert worst_case_separation(list(ACC_COLORS.values())) > worst_case_separation(
        list(ACC_COLORS_LEGACY.values())
    )


def test_assignment_is_the_best_available_from_this_palette():
    """No permutation of the same colours separates better than the one shipped.

    Guards the specific failure that produced this test: the palette was right and
    the class assignment was wrong, which no amount of looking at it revealed.
    """
    shipped = list(ACC_COLORS.values())
    grey = ACC_COLORS["buried_interior"]
    hues = [c for c in shipped if c != grey]
    best = max(
        worst_case_separation(list(order) + [grey])
        for order in itertools.permutations(hues)
    )
    assert worst_case_separation(shipped) >= best - 1e-9


def test_buried_interior_stays_the_least_salient():
    """Grey marks 'nothing to see here' and must not compete with a real class."""
    assert ACC_COLORS["buried_interior"] == "#999999"


def test_palette_remap_is_idempotent_and_total():
    page = "".join(f'<i style="color:{c}"></i>' for c in ACC_COLORS_LEGACY.values())
    once = apply_current_palette(page)
    assert once == apply_current_palette(once)
    for legacy in ACC_COLORS_LEGACY.values():
        assert legacy not in once
    for current in ACC_COLORS.values():
        assert current in once

# --- the portal's own palettes, held to the same bar ------------------------

#: Immune-epitope classes, read from showcase/portal/assets/portal.css. These are
#: drawn side by side in the epitope map, the construct strip and the ribbon, so
#: they must separate from each other exactly as the accessibility classes do.
#: MHC-II was #a98bf5 until 2026-08-14, which measured dE 4.3 against B-cell under
#: protanopia while looking fine (27.1) to normal vision. Nothing tested it.
IMMUNE_CLASSES = {"B-cell": "#5b9df9", "MHC-I": "#4fc47f", "MHC-II": "#8f6cf9"}

#: Pathogen identity. Chrome only — never a data fill — but it shares the screen
#: with both palettes above, so it is held to the same separation.
SPECIES_COLORS = {"pg": "#cbd24b", "tf": "#97d8b1", "td": "#4385b1"}

PORTAL_CSS = (
    pathlib.Path(__file__).resolve().parents[3]
    / "projects" / "YAUVI-PeriodontalPathogens" / "showcase" / "portal" / "assets" / "portal.css"
)


def _token(name: str) -> str:
    """Read a custom property out of the portal stylesheet."""
    for line in PORTAL_CSS.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith(f"--{name}:"):
            return stripped.split(":", 1)[1].strip().rstrip(";")
    raise AssertionError(f"portal.css declares no --{name}")


def test_immune_classes_are_separable():
    """B-cell against MHC-II was the failure this test exists to prevent."""
    separation = worst_case_separation(list(IMMUNE_CLASSES.values()))
    assert separation >= MINIMUM_SEPARATION, (
        f"immune classes collapse to dE {separation:.1f}; B-cell and MHC-II are drawn "
        "side by side in the epitope map, the construct strip and the ribbon."
    )


def test_species_colors_never_collide_with_a_data_palette():
    """Species hue is identity. It must not read as an accessibility or class value."""
    data = list(ACC_COLORS.values()) + list(IMMUNE_CLASSES.values())
    for code, color in SPECIES_COLORS.items():
        for other in data:
            separation = worst_case_separation([color, other])
            assert separation >= MINIMUM_SEPARATION, (
                f"species {code} ({color}) is dE {separation:.1f} from data colour {other}"
            )
    assert worst_case_separation(list(SPECIES_COLORS.values())) >= MINIMUM_SEPARATION


def test_portal_css_matches_the_values_these_tests_derive():
    """Keep the private portal integration in step when that consumer is present."""
    if not PORTAL_CSS.is_file():
        pytest.skip("private portal consumer is outside the standalone reviewer package")
    for name, expected in (("bcell", IMMUNE_CLASSES["B-cell"]),
                           ("mhci", IMMUNE_CLASSES["MHC-I"]),
                           ("mhcii", IMMUNE_CLASSES["MHC-II"]),
                           ("org-pg", SPECIES_COLORS["pg"]),
                           ("org-tf", SPECIES_COLORS["tf"]),
                           ("org-td", SPECIES_COLORS["td"])):
        assert _token(name).lower() == expected.lower(), f"--{name} drifted from the tested value"
