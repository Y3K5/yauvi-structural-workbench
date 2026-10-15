"""The yauvi-fold host bridge that generated pages carry."""
from __future__ import annotations

import pytest

from memorient.viz import (
    FOLD_BRIDGE_PROTOCOL,
    fold_bridge_block,
    fold_bridge_js,
    inject_fold_bridge,
)


PAGE = '<!doctype html><html><body><div id="viewer"></div>' \
       '<script>const v = $3Dmol.createViewer("viewer", {backgroundColor:"white"});</script></body></html>'


def test_bridge_source_declares_the_protocol():
    source = fold_bridge_js()
    assert FOLD_BRIDGE_PROTOCOL in source
    assert "fold-ready" in source and "highlight" in source


def _executable_source() -> str:
    """The bridge with comments stripped, so prose about the design is not tested."""
    source = fold_bridge_js()
    while "/*" in source:
        head, _, rest = source.partition("/*")
        _, _, tail = rest.partition("*/")
        source = head + tail
    return "\n".join(
        line for line in source.splitlines() if not line.strip().startswith("//")
    ).lower()


def test_bridge_carries_no_domain_vocabulary():
    """The page paints what it is told; epitope classes belong to the host.

    Checked against executable source only. The header comment mentions a B-cell
    epitope precisely to say the page does not know what one is, and the legacy
    message name is a compatibility constant rather than domain logic.
    """
    source = _executable_source().replace("oralome-epitope-highlight", "")
    for term in ("b-cell", "mhc-i", "mhc-ii", "epitope_class", "vaccine", "antigen"):
        assert term not in source, f"fold bridge should not know about {term!r}"


def test_the_only_domain_reference_is_the_legacy_compatibility_name():
    """If the legacy alias is ever dropped, this test should be dropped with it."""
    assert "oralome-epitope-highlight" in fold_bridge_js()


def test_injection_adds_the_bridge_before_body_close():
    out = inject_fold_bridge(PAGE)
    assert FOLD_BRIDGE_PROTOCOL in out
    assert out.rstrip().endswith("</body></html>")
    assert out.index(FOLD_BRIDGE_PROTOCOL) < out.rindex("</body>")


def test_injection_normalises_the_viewer_construction():
    """Passing the element id as a string fails on builds expecting a node."""
    out = inject_fold_bridge(PAGE)
    assert '$3Dmol.createViewer(document.getElementById("viewer")' in out
    assert '$3Dmol.createViewer("viewer"' not in out


def test_injection_is_idempotent():
    once = inject_fold_bridge(PAGE)
    assert inject_fold_bridge(once) == once


def test_injection_refuses_a_page_without_a_body():
    with pytest.raises(ValueError, match="no </body>"):
        inject_fold_bridge("<html>no body close</html>")


def test_block_is_self_contained():
    block = fold_bridge_block()
    assert block.count("<script>") == 1 and block.count("</script>") == 1
    assert "yauviFoldLegend" in block


def test_an_older_bridge_is_replaced_not_left_in_place():
    """The reason this function exists: a repair has to reach pages that already
    carry a bridge, or it only ever lands on pages that had none."""
    stale = inject_fold_bridge(PAGE).replace(FOLD_BRIDGE_PROTOCOL, "yauvi-fold/1")
    assert "yauvi-fold/1" in stale

    out = inject_fold_bridge(stale)
    assert FOLD_BRIDGE_PROTOCOL in out
    # One bridge, not two: a second block would leave two message listeners
    # fighting over the same viewer.
    assert out.count("yauviFoldLegend{") == 1
    assert out.count('id = "yauviFoldLegend"') == 1
    assert out == inject_fold_bridge(PAGE)


def test_a_version_marker_with_no_recognisable_block_is_refused():
    """Better to stop than to append a second listener to an unknown layout."""
    with pytest.raises(ValueError, match="could not be located"):
        inject_fold_bridge(
            '<!doctype html><html><body><script>var p = "yauvi-fold/1";</script></body></html>'
        )


def test_older_protocol_messages_are_still_accepted():
    """Pages regenerate independently of the hosts driving them, and one campaign's
    workbench still speaks /1."""
    source = fold_bridge_js()
    assert '"yauvi-fold/1": true' in source and '"yauvi-fold/2": true' in source


def test_highlights_are_spheres_over_an_untouched_cartoon():
    """Hue is spent on accessibility. A highlight that recoloured the cartoon
    underneath spent that channel a second time."""
    source = _executable_source()
    # addStyle keeps the cartoon; setStyle would replace it.
    assert "addstyle" in source
    assert "cartoon: { color: color }" not in source


def test_a_highlight_marks_the_span_even_without_a_surface():
    """Surface generation is the one call here that can fail on an old build, and a
    highlight that silently renders nothing is worse than a plain one."""
    source = _executable_source()
    stick = source.index("view.addstyle(selection, { stick:")
    # the call site, not the definition above it
    patch = source.index("addpatch(selection, colour);")
    assert stick < patch, "the fallback sticks are drawn after the surface, not before"
    assert "typeof view.addsurface !== \"function\"" in source


def test_surfaces_are_removed_by_handle():
    """A surface is an object on the viewer, not a style on the atoms, so restyling
    the residues does not clear it."""
    source = _executable_source()
    assert "removesurface" in source and "droppatches" in source
    assert source.index("droppatches();") < source.index("painted.concat(inspected)")


def test_the_camera_centres_rather_than_zooming_to_a_highlight():
    """zoomTo frames the span alone, which fills the viewport with empty space and
    pushes the protein out of sight. The highlight has to stay on something."""
    source = _executable_source()
    assert "view.center(selection" in source
    assert "view.zoomto(selection" not in source
