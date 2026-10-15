/* yauvi-fold/4 — the postMessage protocol between a generated fold page and its host.
 *
 * A fold page is a self-contained 3Dmol.js document in an iframe. The host needs
 * to highlight residue ranges on it, be told when it is ready, and follow the
 * pointer for a linked sequence track. That contract used to live as a string
 * inside one campaign's build script and was patched into generated pages after
 * the fact, so the page and the host could drift apart silently.
 *
 * This is the single copy, kept beside membrane_bilayer.js for the same reason.
 *
 * The protocol is deliberately free of any vaccine vocabulary. This page knows
 * about residues, colours and camera state; it does not know what a "B-cell
 * epitope" is. The host sends the colour it wants and the label to display, so
 * the same bridge serves any consumer.
 *
 *   page   -> host : {protocol, type:"fold-ready", path}
 *                    {protocol, type:"residue-hover", resid, chain}
 *   host   -> page : {protocol, type:"highlight", ranges:[...], inspected:[...], caption}
 *                    {protocol, type:"view-sync", view}
 *
 * A range is {start, end, color, labelColor?, label?, chain?, style?}.
 *
 * WHAT CHANGED IN /4 — a discrete range is drawn as a translucent patch over the
 * residues plus sticks beneath it, and the camera turns to face it.
 *
 * Beads floated above the fold: legible, but they read as objects added to the
 * molecule rather than as a region of it. A patch reads as part of the surface,
 * which is what an epitope is. The colour moved off white for the same reason it
 * moved off class hue — white competes with the pale end of the accessibility
 * palette, and a saturated red sits clear of every one of the six classes,
 * including the vermillion used for periplasmic.
 *
 * Sticks are drawn first and unconditionally. Surface generation is the one call
 * here that can be slow or fail on an old build, and a highlight that silently
 * renders nothing is worse than a plain one.
 *
 * `focus: false` on a range suppresses the camera move, for a host repainting
 * without the reader having asked to go anywhere.
 *
 * WHAT CHANGED IN /3 — a range may ask for `style: "surface"` instead of the
 * default bead.
 *
 * A bead per residue is the right mark for a peptide: a dozen or so residues,
 * read as a path along the backbone. It is the wrong mark for a region. A host
 * showing "everything on this protein that is antibody-accessible" is sending
 * hundreds of residues in dozens of runs, and that many spheres is an object that
 * hides the fold rather than a highlight on it. `surface` draws those as sticks,
 * which thickens the existing cartoon instead of covering it, so a region and a
 * peptide stay visually distinct even when drawn in the same colour.
 *
 * WHAT CHANGED IN /2 — highlighted ranges are drawn as spheres on CA atoms over an
 * untouched cartoon, not as recoloured cartoon plus thick sticks.
 *
 * Two reasons. First, the generated page already spends hue on accessibility, one
 * class per colour; recolouring the cartoon underneath a highlight put a second
 * meaning on the same channel, so a highlighted residue stopped reporting whether
 * it was antibody-accessible. Second, a stick overlay reads as a thicker part of
 * the same object. One bead per residue reads as a separate object laid along the
 * backbone, which is what tracing a peptide onto a barrel actually needs, and it
 * stays legible when several ranges are lit at once.
 *
 * The page still paints whatever colour the host asks for. A host wanting the
 * beads to carry meaning sends a colour per range; a host wanting them uniform
 * sends the same one, and passes `labelColor` so the legend can keep a hue the
 * geometry no longer carries.
 *
 * `inspected` ranges are the residue under the pointer: looked at, not chosen.
 * They are drawn as a larger bead rather than a different colour, because white
 * is now the selected channel for at least one host. Restoring a range returns it
 * to the page's own per-residue colouring, which is this page's data, not the
 * host's.
 *
 * /1 messages are still accepted. Pages are regenerated independently of the hosts
 * that drive them, and one campaign's workbench still speaks /1.
 */
(function () {
  "use strict";

  var PROTOCOL = "yauvi-fold/4";
  /* Accepted so a host that has not been updated keeps working. */
  var ACCEPTED = {
    "yauvi-fold/4": true, "yauvi-fold/3": true, "yauvi-fold/2": true, "yauvi-fold/1": true
  };
  var LEGACY_HIGHLIGHT = "oralome-epitope-highlight";

  /* Bead radii in ångström. CA-CA spacing is ~3.8 Å, so SELECTED leaves a visible
     gap between consecutive residues — a chain of beads rather than a tube — and
     INSPECTED overlaps its neighbours enough to read as one residue picked out. */
  var SELECTED_RADIUS = 1.6;
  var INSPECTED_RADIUS = 2.3;
  /* Thick enough to find against the cartoon, thin enough that a few hundred of
     them still read as the same molecule. */
  var SURFACE_RADIUS = 0.5;
  /* Sits clear of all six accessibility classes, vermillion included. */
  var MARK = "#ff2d55";
  var MARK_STICK_RADIUS = 0.42;
  var PATCH_OPACITY = 0.72;
  var FOCUS_MS = 650;

  /* Surfaces are objects on the viewer, not styles on atoms, so they have to be
     removed by handle rather than cleared by restyling the residues. */
  var patches = [];

  function dropPatches() {
    for (var i = 0; i < patches.length; i += 1) {
      try { view.removeSurface(patches[i]); } catch (error) { /* already gone */ }
    }
    patches = [];
  }

  function addPatch(selection, colour) {
    if (typeof view.addSurface !== "function") { return; }
    try {
      var made = view.addSurface(
        ($3Dmol && $3Dmol.SurfaceType && $3Dmol.SurfaceType.VDW) || 1,
        { opacity: PATCH_OPACITY, color: colour },
        selection
      );
      /* 2.x returns a promise of the handle; 1.x returns the handle itself. */
      if (made && typeof made.then === "function") {
        made.then(function (handle) { patches.push(handle); view.render(); });
      } else if (made !== undefined) {
        patches.push(made);
      }
    } catch (error) { /* older build: the sticks below still mark the span */ }
  }

  /* `v` and `colors` are top-level const bindings, so they are intentionally not
     properties of window. Resolve the lexical bindings first; the old window-only
     lookup selected the #viewer element and made epitope tracing fail at render(). */
  var view = typeof v !== "undefined" && v && typeof v.render === "function" ? v : null;
  var ledger = typeof colors !== "undefined" && Array.isArray(colors) ? colors : [];
  if (!view) { return; }

  var byResidue = {};
  for (var i = 0; i < ledger.length; i += 1) { byResidue[ledger[i].resid] = ledger[i]; }

  var painted = [];
  var inspected = [];

  function send(message) {
    message.protocol = PROTOCOL;
    try { window.parent.postMessage(message, "*"); } catch (error) { /* not framed */ }
  }

  /* setStyle, not addStyle: it replaces every style on the matching atoms, which is
     what removes the bead. The selection covers the whole residue, so the CA the
     bead sat on is included. */
  function restore(range) {
    var chain = range.chain || "A";
    for (var resid = range.start; resid <= range.end; resid += 1) {
      var row = byResidue[resid];
      if (!row) {
        /* No per-residue ledger on this page: clear back to the base style. */
        view.setStyle({ chain: chain, resi: resid }, { cartoon: {} });
        continue;
      }
      view.setStyle({ chain: row.chain, resi: row.resid }, { cartoon: { color: row.color } });
      if (row.extracellular) {
        view.addStyle({ chain: row.chain, resi: row.resid },
                      { stick: { color: row.color, radius: 0.2 } });
      }
    }
  }

  /* addStyle, so the cartoon and its accessibility colour survive underneath. */
  function paint(range) {
    var selection = { chain: range.chain || "A", resi: range.start + "-" + range.end };
    if (range.style === "surface") {
      /* No atom filter: the whole residue thickens, which is what makes a region
         read as a region rather than as a very long peptide. */
      view.addStyle(selection,
        { stick: { color: range.color || "#ffffff", radius: SURFACE_RADIUS } });
      return;
    }
    if (range.style === "bead") {
      selection.atom = "CA";
      view.addStyle(selection,
        { sphere: { color: range.color || MARK, radius: SELECTED_RADIUS } });
      return;
    }
    /* Default: a patch over the span, with sticks under it so the mark survives a
       viewer that cannot build a surface. */
    var colour = range.color || MARK;
    view.addStyle(selection, { stick: { color: colour, radius: MARK_STICK_RADIUS } });
    addPatch(selection, colour);
    if (range.focus !== false && typeof view.center === "function") {
      /* Centre on it, do not zoom to it. zoomTo frames the selection alone, and a
         fifteen-residue span framed alone fills the viewport with empty space and
         pushes the rest of the protein out of sight — the reader loses the object
         the highlight was supposed to be *on*. Centring keeps the whole fold and
         moves the marked region to the middle. */
      try { view.center(selection, FOCUS_MS); } catch (error) { /* older build */ }
    }
  }

  function paintInspected(range) {
    view.addStyle(
      { chain: range.chain || "A", resi: range.start + "-" + range.end, atom: "CA" },
      { sphere: { color: range.color || "#ffffff", radius: INSPECTED_RADIUS } }
    );
  }

  function legendHost() {
    var host = document.getElementById("yauviFoldLegend");
    if (host) { return host; }
    host = document.createElement("div");
    host.id = "yauviFoldLegend";
    host.setAttribute("aria-live", "polite");
    document.body.appendChild(host);
    return host;
  }

  /* The caption is whatever the host sent. This page does not invent wording for
     a domain it knows nothing about.

     Only discrete ranges are listed. A `surface` range is one run of a region that
     may have a hundred runs in it, and enumerating those turned the legend into a
     panel that covered the structure it was annotating — the count belongs in the
     caption the host already sends, not in a hundred labels. Even the discrete
     list is capped: past a few, a legend stops being a key and becomes a table. */
  var LEGEND_MAX = 4;

  function describe(caption) {
    var host = legendHost();
    if (!painted.length && !inspected.length) {
      host.innerHTML = "<strong>" + (caption || "Nothing highlighted") + "</strong>";
      return;
    }
    var discrete = [];
    for (var i = 0; i < painted.length; i += 1) {
      if (painted[i].style !== "surface") { discrete.push(painted[i]); }
    }
    var parts = discrete.slice(0, LEGEND_MAX).map(function (range) {
      var label = range.label || (range.start + "-" + range.end);
      /* labelColor lets a host draw uniform beads and still colour the legend. */
      var swatch = range.labelColor || range.color || "#fff";
      return '<span style="color:' + swatch + '">' + label + "</span>";
    });
    if (discrete.length > LEGEND_MAX) {
      parts.push('<span class="more">+' + (discrete.length - LEGEND_MAX) + " more</span>");
    }
    if (inspected.length) {
      parts.push('<span class="inspect">inspecting ' + inspected[0].start + "-" +
                 inspected[0].end + " · not chosen</span>");
    }
    host.innerHTML = "<strong>" + (caption || (painted.length + " highlighted")) + "</strong>" +
                     parts.join("");
  }

  function applyHighlight(data) {
    dropPatches();
    painted.concat(inspected).forEach(restore);
    painted = Array.isArray(data.ranges) ? data.ranges : [];
    inspected = Array.isArray(data.inspected) ? data.inspected : [];
    painted.forEach(paint);
    inspected.forEach(paintInspected);
    describe(data.caption);
    view.render();
  }

  window.addEventListener("message", function (event) {
    if (event.source !== window.parent || !event.data) { return; }
    var data = event.data;
    var isCurrent = ACCEPTED[data.protocol] === true;
    var isLegacy = data.type === LEGACY_HIGHLIGHT;
    if (!isCurrent && !isLegacy) { return; }

    if (isLegacy || data.type === "highlight") { applyHighlight(data); return; }
    if (data.type === "view-sync" && data.view) {
      view.setView(data.view);
      view.render();
    }
  });

  /* Camera persistence, so switching targets and coming back does not reset the
     view the reader had set up. */
  var viewKey = "yauvi-fold-view:" + location.pathname;
  try {
    var stored = JSON.parse(sessionStorage.getItem(viewKey) || "null");
    if (stored) { view.setView(stored); view.render(); document.body.dataset.viewState = "restored"; }
  } catch (error) { /* storage unavailable */ }

  function saveView() {
    try {
      sessionStorage.setItem(viewKey, JSON.stringify(view.getView()));
      document.body.dataset.viewState = "saved";
    } catch (error) { /* storage unavailable */ }
  }
  document.addEventListener("pointerup", saveView);
  document.addEventListener("wheel", saveView, { passive: true });

  /* Pointer -> residue, for a host that draws a linked sequence track. */
  if (typeof view.setHoverable === "function") {
    try {
      view.setHoverable({}, true, function (atom) {
        if (atom && atom.resi !== undefined) {
          send({ type: "residue-hover", resid: atom.resi, chain: atom.chain || "A" });
        }
      }, function () {
        send({ type: "residue-hover", resid: null, chain: null });
      });
    } catch (error) { /* older 3Dmol build */ }
  }

  send({ type: "fold-ready", path: location.pathname });
}());
