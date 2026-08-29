"""memorient.viz — 3D export for oriented structures.

Two exporters, both driven by an :class:`~memorient.orientor.OrientationResult`:

* :func:`display_oriented` → a JSON-serializable descriptor a 3Dmol.js viewer can consume:
  the oriented PDB text, a per-residue colour map keyed by accessibility, and a **membrane
  slab descriptor** (leaflet bounds + thickness) that is present only when the context has a
  bilayer — a soluble protein gets no slab, matching the "metrics are not a blanket" rule.
* :func:`write_pymol_script` → a ``.pml`` that loads the oriented PDB, colours residues by
  membrane zone / accessibility, and draws the two leaflet planes as pseudo-atoms.

The oriented frame produced by the orientor already places the membrane centre at the origin
with the extracellular side at +Z, so the viewer just needs the leaflet z-bounds.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from .geometry import to_pdb_string
from .membrane import (
    ACC_ANTIBODY,
    ACC_BURIED,
    ACC_LIPID,
    ACC_LPS_SHIELDED,
    ACC_PERIPLASMIC,
    ACC_PORE,
    LPS_BUFFER,
)

# Colour scheme shared by both exporters (hex, viewer-agnostic).
#
# Okabe-Ito, with the class assignment chosen by search rather than by eye. The
# selection criterion was the worst-case CIE76 dE between any two classes under
# normal, deuteranope, protanope and tritanope vision (Vienot-Brettel-Mollon
# 1999 simulation); `tests/test_palette.py` re-derives it and fails if the
# assignment stops being optimal.
#
# Measured worst-case separation:
#     previous scheme            13.7
#     Okabe-Ito, assigned by eye  7.6   (reddish purple collapsed toward grey)
#     this assignment            14.9
#
# The by-eye attempt was worse than what it replaced, which is why this is
# measured. Note the previous scheme was not catastrophic on this metric; its
# real weakness was that its three most consequential classes sat along the
# red/orange/yellow axis, so the failure mode was concentrated exactly where the
# reader's decision is made.
#
# Hue carries accessibility class and nothing else. Epitope highlighting uses
# geometry — one sphere per CA over an untouched cartoon — rather than more hues,
# so the two meanings never compete for the same channel. An earlier bridge
# recoloured the cartoon under a highlight, which silently spent this channel
# twice; see fold_bridge.js.
ACC_COLORS: Dict[str, str] = {
    ACC_ANTIBODY: "#E69F00",       # orange — the publishable epitope surface
    ACC_LPS_SHIELDED: "#56B4E9",   # sky blue — extracellular but LPS-buried
    ACC_LIPID: "#F0E442",          # yellow — lipid-facing
    ACC_PORE: "#0072B2",           # blue — lumen/water-facing
    ACC_PERIPLASMIC: "#D55E00",    # vermillion — periplasmic
    ACC_BURIED: "#999999",         # grey — buried interior, deliberately low salience
}

#: The scheme used before the Okabe-Ito change. Kept so pages generated earlier
#: stay readable and can be remapped by :func:`apply_current_palette`, and so an
#: archived figure's colours can still be interpreted.
ACC_COLORS_LEGACY: Dict[str, str] = {
    ACC_ANTIBODY: "#e6194b",
    ACC_LPS_SHIELDED: "#f58231",
    ACC_LIPID: "#ffe119",
    ACC_PORE: "#4363d8",
    ACC_PERIPLASMIC: "#911eb4",
    ACC_BURIED: "#a9a9a9",
}

ZONE_COLORS: Dict[str, str] = {
    "extracellular": "#0072B2",
    "extracellular_interface": "#56B4E9",
    "hydrophobic_core": "#E69F00",
    "periplasmic_interface": "#009E73",
    "periplasmic": "#CC79A7",
}

ZONE_COLORS_LEGACY: Dict[str, str] = {
    "extracellular": "#e6194b",
    "extracellular_interface": "#f58231",
    "hydrophobic_core": "#ffe119",
    "periplasmic_interface": "#42d4f4",
    "periplasmic": "#911eb4",
}


def apply_current_palette(html: str) -> str:
    """Rewrite a generated page's legacy accessibility hexes to the current scheme.

    Orientation results are expensive to recompute, so pages generated under the
    old palette are remapped in place rather than regenerated. The mapping is
    exact-hex and class-wise, so a page that already uses the current scheme is
    unchanged and the operation is idempotent.
    """
    for key, legacy in list(ACC_COLORS_LEGACY.items()) :
        current = ACC_COLORS[key]
        html = html.replace(legacy, current).replace(legacy.upper(), current)
    for key, legacy in ZONE_COLORS_LEGACY.items():
        html = html.replace(legacy, ZONE_COLORS[key]).replace(legacy.upper(), ZONE_COLORS[key])
    return html


#: The generated template loads 3Dmol from the project CDN, which is correct for a
#: standalone page a user opens with a network. It is wrong for anything shipped:
#: over ``file://``, inside the ``.app`` bundle, or on a machine with no route out,
#: the script never arrives and the viewer renders an empty div — with no error a
#: reader would notice. Hosts that ship these pages call `vendor_3dmol()`.
_CDN_3DMOL = "https://3Dmol.org/build/3Dmol-min.js"


def vendor_3dmol(html: str, relative_path: str) -> str:
    """Point a generated page at a locally vendored 3Dmol instead of the CDN.

    Idempotent: a page already vendored is returned unchanged.
    """
    return html.replace(_CDN_3DMOL, relative_path)


#: Light chrome baked into `_HTML_TEMPLATE`, and the dark value each maps to. The
#: page was designed to be opened on its own, where white is right; embedded in a
#: dark host it reads as a lit panel punched through the page.
_DARK_CHROME = (
    ("background:#fafafa", "background:#0a0f14"),
    ('backgroundColor:"white"', 'backgroundColor:"#0a0f14"'),
    ("border-bottom:1px solid #ddd", "border-bottom:1px solid rgba(168,188,199,0.14)"),
    ("#head span{color:#666", "#head span{color:#74899a"),
    # The unassigned cartoon. Light grey on a dark ground out-shouts every
    # accessibility hue on top of it, which inverts the figure's whole point.
    ('cartoon:{color:"#dddddd"}', 'cartoon:{color:"#33424f"}'),
)


def apply_dark_chrome(html: str) -> str:
    """Recolour a generated page for a dark host, leaving the data colours alone.

    Only chrome moves: background, rule, caption and the *unassigned* cartoon. Every
    per-residue accessibility colour is untouched, so the figure still says exactly
    what it said before. Idempotent, like `apply_current_palette`.
    """
    for light, dark in _DARK_CHROME:
        html = html.replace(light, dark)
    if "color:#e8f0f4" not in html:
        html = html.replace(
            "body{font-family:", "body{color:#e8f0f4;font-family:", 1
        )
    return html


def _membrane_slab(result) -> Optional[Dict[str, object]]:
    """Leaflet bounds for the oriented frame, or ``None`` when the context has no bilayer."""
    from .contexts import get_context

    ctx = get_context(result.context)
    if not ctx.has_bilayer or result.fit is None:
        return None
    d = float(result.fit.half_thickness)
    slab = {
        "half_thickness": d,
        "core_lower_z": -d,
        "core_upper_z": d,
        "extracellular_z": d,       # +Z is extracellular (orientor convention)
        "periplasmic_z": -d,
        "asymmetric": ctx.is_asymmetric,
    }
    if ctx.is_asymmetric:
        # LPS occupies the outer leaflet + a buffer band above the extracellular interface
        slab["lps_upper_z"] = d + LPS_BUFFER
        slab["lps_leaflet"] = "extracellular"
    return slab


def display_oriented(result) -> Dict[str, object]:
    """Build a 3Dmol.js-ready descriptor from an :class:`OrientationResult`."""
    labels = result.labels.labels
    residue_colors: List[Dict[str, object]] = []
    for l in labels:
        color = ACC_COLORS.get(l.accessibility, "#a9a9a9")
        residue_colors.append({
            "resid": l.resid,
            "chain": l.chain,
            "accessibility": l.accessibility,
            "zone": l.zone,
            "extracellular": bool(l.extracellular),
            "color": color,
        })
    return {
        "context": result.context,
        "method": result.method,
        "label": result.label,
        "pdb": to_pdb_string(result.structure),
        "orientation": {"extracellular_axis": "+z", "membrane_center": [0.0, 0.0, 0.0]},
        "membrane_slab": _membrane_slab(result),   # None for soluble/anchored
        "residue_colors": residue_colors,
        "surface_set": list(result.labels.surface_set),
        "color_legend": ACC_COLORS,
    }


def write_pymol_script(result, path: str) -> None:
    """Write a PyMOL ``.pml`` that colours the oriented structure by accessibility + draws the slab."""
    lines: List[str] = []
    lines.append("# memorient oriented view")
    lines.append(f"# context={result.context} method={result.method} label={result.label}")
    lines.append("load oriented.pdb, mol")
    lines.append("hide everything, mol")
    lines.append("show cartoon, mol")
    lines.append("bg_color white")
    lines.append("set cartoon_transparency, 0.1")

    # group residues by accessibility and colour each group
    by_acc: Dict[str, List[int]] = {}
    for l in result.labels.labels:
        by_acc.setdefault(l.accessibility, []).append(l.resid)
    for acc, resids in by_acc.items():
        hexc = ACC_COLORS.get(acc, "#a9a9a9")
        r, g, b = (int(hexc[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
        cname = f"acc_{acc}"
        sel = "+".join(str(x) for x in resids)
        lines.append(f"set_color {cname}, [{r:.3f}, {g:.3f}, {b:.3f}]")
        lines.append(f"select {acc}, mol and resi {sel}")
        lines.append(f"color {cname}, {acc}")

    # highlight the antibody-accessible epitope surface as sticks
    if result.labels.surface_set:
        sel = "+".join(str(x) for x in result.labels.surface_set)
        lines.append(f"select epitope_surface, mol and resi {sel}")
        lines.append("show sticks, epitope_surface")

    # membrane slab as two translucent planes (pseudo-atoms), only when there is a bilayer
    slab = _membrane_slab(result)
    if slab is not None:
        d = slab["half_thickness"]
        lines.append(f"pseudoatom mem_ec, pos=[0,0,{d:.1f}]")
        lines.append(f"pseudoatom mem_peri, pos=[0,0,{-d:.1f}]")
        lines.append("# membrane core spans z = [%.1f, %.1f]; extracellular is +z" % (-d, d))
    lines.append("orient mol")
    lines.append("zoom mol, 5")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


_HTML_TEMPLATE = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>memorient — {title}</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
 body{{font-family:-apple-system,Helvetica,Arial,sans-serif;margin:0;background:#fafafa}}
 #head{{padding:10px 16px;border-bottom:1px solid #ddd}}
 #head b{{font-size:15px}} #head span{{color:#666;font-size:13px}}
 #viewer{{width:100%;height:78vh;position:relative}}
 #legend{{padding:8px 16px;font-size:12px}} .sw{{display:inline-block;width:12px;height:12px;
   margin:0 4px -1px 12px;border-radius:2px}}
</style></head><body>
<div id="head"><b>{title}</b> &nbsp; <span>{subtitle}</span></div>
<div id="viewer"></div>
<div id="legend">{legend}</div>
<script>{bilayer_js}</script>
<script>
const pdb = {pdb_json};
const colors = {colors_json};
const slab = {slab_json};
const v = $3Dmol.createViewer("viewer", {{backgroundColor:"white"}});
v.addModel(pdb, "pdb");
v.setStyle({{}}, {{cartoon:{{color:"#dddddd"}}}});
for (const r of colors) {{
  v.setStyle({{resi:r.resid}}, {{cartoon:{{color:r.color}}}});
  if (r.extracellular) v.addStyle({{resi:r.resid}}, {{stick:{{color:r.color,radius:0.2}}}});
}}
// The bilayer is drawn by the shared MembraneBilayer module (inlined above), so this
// page and the viewers that vendor this module cannot drift apart.
// +Z is extracellular by the orientor convention, so the leaflets are z-planes.
if (slab) {{
  MembraneBilayer.draw(v, {{
    axis: "z",
    core: [slab.core_lower_z, slab.core_upper_z],
    leaflets: [
      {{at: slab.core_upper_z, color: "#6d5bd0"}},   // extracellular core boundary
      {{at: slab.core_lower_z, color: "#6d5bd0"}},   // periplasmic core boundary
      // LPS keeps its own hue: in a gram-negative outer membrane the outer leaflet is
      // not the same chemistry as the inner one, and one colour would hide that.
      slab.lps_upper_z ? {{at: slab.lps_upper_z, color: "#c8791f"}} : null
    ].filter(Boolean)
    // No haze: this page renders on WHITE, where even a 7%-opacity core slab reads as a
    // pale rectangle with hard edges. The leaflet planes already bound the core.
  }});
}}
v.zoomTo();
// 3Dmol's default camera looks straight down -Z, which is the membrane normal here — so
// the bilayer would be seen face-on as concentric rings. Tilt to a side view so it reads
// as a membrane, with extracellular (+Z) up and periplasm down.
// -90, not +90: a +90 x-rotation puts +Z at the BOTTOM, which would show the LPS leaflet
// and the antibody-accessible residues below the periplasmic domain — an inverted
// membrane in an orientation tool is worse than no membrane at all.
if (slab) v.rotate(-90, "x");
v.render();
</script></body></html>
"""


def _bilayer_js() -> str:
    """The shared membrane module's source, for inlining into a generated page.

    Read from disk rather than embedded as a Python string so there is exactly one copy to
    edit. Viewers outside this distribution vendor a byte-identical copy; this file is the
    canonical source for all of them.
    """
    from pathlib import Path
    return (Path(__file__).with_name("membrane_bilayer.js")).read_text(encoding="utf-8")


FOLD_BRIDGE_PROTOCOL = "yauvi-fold/4"

#: Matches an already-injected bridge block — the legend styling through the end of
#: the script that follows it. Used to replace an older bridge rather than leave it
#: in place; see :func:`inject_fold_bridge`.
_FOLD_BRIDGE_BLOCK_RE = re.compile(
    r"\n?<style>\s*#yauviFoldLegend.*?</style>\s*<script>.*?</script>\s*",
    re.DOTALL,
)

#: Any bridge version marker, so an older one can be recognised and replaced.
_FOLD_BRIDGE_VERSION_RE = re.compile(r"yauvi-fold/(\d+)")

#: Minimal styling for the legend the bridge maintains. Kept here rather than in
#: the bridge so a host page can restyle it without forking the protocol.
_FOLD_BRIDGE_CSS = """
<style>
#yauviFoldLegend{position:fixed;right:10px;top:10px;z-index:20;max-width:42%;max-height:30%;
overflow:hidden;padding:7px 9px;
background:rgba(7,17,27,.90);border:1px solid rgba(255,255,255,.22);border-radius:7px;color:#edf5ff;
font:11px/1.35 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;pointer-events:none}
#yauviFoldLegend strong{display:block;margin-bottom:3px;color:#fff}
#yauviFoldLegend span{display:inline-block;margin-right:8px}
#yauviFoldLegend .inspect{color:#fff}
#yauviFoldLegend .more{opacity:.7}
</style>
"""


def fold_bridge_js() -> str:
    """Source of the ``yauvi-fold/1`` host bridge, for inlining into a page.

    Read from disk for the same reason as :func:`_bilayer_js`: exactly one copy to
    edit. Exposed publicly because a builder may need to add the bridge to a page
    it generated itself, or to one generated before the bridge existed.
    """
    from pathlib import Path
    return (Path(__file__).with_name("fold_bridge.js")).read_text(encoding="utf-8")


def fold_bridge_block() -> str:
    """The bridge as a drop-in HTML fragment: legend styling plus the script."""
    return f"{_FOLD_BRIDGE_CSS}<script>\n{fold_bridge_js()}\n</script>\n"


def inject_fold_bridge(html: str) -> str:
    """Return ``html`` with the current fold bridge in it, replacing an older one.

    Also normalises the 3Dmol viewer construction: passing the element id as a
    string fails on builds that expect a node, and every consumer was patching
    that separately.

    This used to return early whenever *any* bridge was present. That made the
    injection idempotent, which was the intent, but it also made it inert: every
    fix to the bridge stopped at the pages that had none, and the generated pages
    diverged from this file without anything reporting it. Version the block
    instead — same page in, same page out when it is already current, replaced
    when it is behind.
    """
    html = html.replace(
        '$3Dmol.createViewer("viewer",',
        '$3Dmol.createViewer(document.getElementById("viewer"),',
    )
    found = _FOLD_BRIDGE_VERSION_RE.search(html)
    if found:
        current = int(_FOLD_BRIDGE_VERSION_RE.search(FOLD_BRIDGE_PROTOCOL).group(1))
        if int(found.group(1)) >= current:
            return html
        # The block already opens with its own newline, and the pattern consumed the
        # one that preceded it, so an upgraded page is byte-identical to a fresh one.
        replaced, count = _FOLD_BRIDGE_BLOCK_RE.subn(
            lambda _: fold_bridge_block(), html, count=1
        )
        if count != 1:
            # A page carrying a version marker but no recognisable block is not
            # something to guess at: appending a second bridge would leave two
            # message listeners fighting over the same viewer.
            raise ValueError(
                f"page declares {found.group(0)} but its bridge block could not be "
                "located, so it cannot be upgraded in place"
            )
        return replaced
    if "</body>" not in html:
        raise ValueError("cannot inject the fold bridge: page has no </body>")
    head, _, tail = html.rpartition("</body>")
    return head + fold_bridge_block() + "</body>" + tail


def write_3dmol_html(result, path: str) -> None:
    """Write a self-contained 3Dmol.js HTML page showing the oriented structure + membrane slab.

    Cartoon is coloured by accessibility (same palette as the PyMOL export); extracellular
    residues are additionally drawn as sticks. The membrane is drawn by the shared
    MembraneBilayer module, inlined into the page so it stays openable over file:// with no
    sibling fetches. +Z is extracellular by the orientor convention, so top = outside.
    """
    import json as _json

    disp = display_oriented(result)
    legend_bits = []
    for acc, hexc in ACC_COLORS.items():
        if any(r["accessibility"] == acc for r in disp["residue_colors"]):
            legend_bits.append(f'<span class="sw" style="background:{hexc}"></span>{acc}')
    sm = result.summary()
    subtitle = (f"context={result.context} · method={result.method} · {result.label} · "
                f"conf={sm.get('confidence')} · {sm.get('n_extracellular', 0)} extracellular · "
                f"surface set {len(result.labels.surface_set)}")
    html = _HTML_TEMPLATE.format(
        title=f"{result.label} — oriented",
        subtitle=subtitle,
        legend="".join(legend_bits),
        pdb_json=_json.dumps(disp["pdb"]),
        colors_json=_json.dumps(disp["residue_colors"]),
        slab_json=_json.dumps(disp["membrane_slab"]),
        bilayer_js=_bilayer_js(),
    )
    # Newly generated pages carry the bridge, so nothing has to patch them later.
    html = inject_fold_bridge(html)
    with open(path, "w") as fh:
        fh.write(html)
