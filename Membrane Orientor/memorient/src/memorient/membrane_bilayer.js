/* CANONICAL SOURCE — membrane-bilayer.js
 *
 * Home: memorient, because this is how the membrane orientor draws a membrane.
 * Every viewer we own inlines THIS file at build time; none of them link to it at
 * runtime, because each generated page has to stay self-contained and file://-safe.
 *
 * Consumers (kept byte-identical by sync_membrane_bilayer.py --check):
 *   - memorient           src/memorient/viz.py            -> oriented.html
 *   - TDVax-YAUVI portal  yauvi/assets/membrane-bilayer.js -> Targets view
 *   - Triple Vax portal   src/redvax/viz/build_membrane_portal.py -> Membrane Atlas
 *
 * Do not edit a vendored copy. Edit this file, then run:
 *   python3 sync_membrane_bilayer.py
 */

/* membrane-bilayer.js — one way to draw a membrane, for every 3Dmol viewer we own.
 *
 * The bilayer is drawn the way a structure viewer shows one (OPM / AlphaFold-DB style
 * dummy head groups): leaflets of small beads rather than solid slabs, with the
 * protein's own footprint cut out, so the barrel stays readable THROUGH the membrane
 * instead of sitting behind a tinted box.
 *
 * Written as a module because the same drawing exists in three viewers — this portal,
 * memorient's oriented.html generator, and the Triple Vax membrane portal — and they had
 * begun to drift. It is deliberately free of any YAUVI assumption:
 *
 *   - `axis` picks the membrane normal, so a Y-normal viewer (Triple Vax) and a Z-normal
 *     viewer (this portal, memorient) share one implementation.
 *   - `leaflets` is a list, so a symmetric bilayer, an asymmetric OM with a separate LPS
 *     band, or a single anchor plane are all the same call with different data.
 *   - `haze` is optional: it reads well on a dark background and shows as a pale
 *     rectangle with hard edges on white, so the caller decides.
 *
 * No dependencies, no build step, no module system — a plain global, matching how the
 * rest of this portal's assets load.
 *
 * Usage:
 *   MembraneBilayer.draw(viewer, {
 *     axis: "z",
 *     core: [coreMinZ, coreMaxZ],
 *     leaflets: [{ at: outerZ, color: "#f59e0b" }, { at: innerZ, color: "#06b6d4" }],
 *     haze: { color: "#64748b", opacity: 0.07 },
 *     selection: { chain: "A" }
 *   });
 *   // -> { cx, cy, radius, beads }   (cx/cy are in the two in-plane axes)
 */
(function (root) {
  "use strict";

  // Which coordinate keys are in-plane (u, v) and which is the normal (n).
  var AXES = {
    z: { u: "x", v: "y", n: "z" },
    y: { u: "x", v: "z", n: "y" },
    x: { u: "y", v: "z", n: "x" }
  };

  var DEFAULTS = {
    axis: "z",
    clearance: 1.2,     // gap between the protein's radius and the first bead; null = no hole
    extent: 18,         // how far past the protein the bead field reaches
    beads: 130,         // target bead count per leaflet — step is solved from this
    radius: 1.5,        // bead radius
    opacity: 0.85,
    jitter: 0.3,        // in-plane jitter as a fraction of the lattice step
    thickness: 1.6,     // out-of-plane jitter, so a leaflet is not a perfect sheet
    selection: {}
  };

  function opt(spec, key) {
    return spec[key] === undefined || spec[key] === null && key !== "clearance"
      ? DEFAULTS[key]
      : spec[key];
  }

  /* Footprint of the protein IN THE MEMBRANE PLANE.
   *
   * Measured only from atoms between the core bounds when those are given. Using the whole
   * protein sizes the hole to whatever hangs outside the bilayer — BamA's periplasmic POTRA
   * arm, for one — and blows the membrane up to several times the barrel it belongs to. */
  function footprint(viewer, spec) {
    var ax = AXES[opt(spec, "axis")] || AXES.z;
    var model = typeof viewer.getModel === "function" ? viewer.getModel() : null;
    var atoms = model && typeof model.selectedAtoms === "function"
      ? model.selectedAtoms(opt(spec, "selection"))
      : [];
    if (!atoms.length) return null;

    var core = spec.core;
    var inCore = core
      ? atoms.filter(function (a) { return a[ax.n] >= core[0] && a[ax.n] <= core[1]; })
      : atoms;
    if (!inCore.length) inCore = atoms;

    var us = inCore.map(function (a) { return a[ax.u]; });
    var vs = inCore.map(function (a) { return a[ax.v]; });
    var cx = (Math.min.apply(Math, us) + Math.max.apply(Math, us)) / 2;
    var cy = (Math.min.apply(Math, vs) + Math.max.apply(Math, vs)) / 2;
    var r = 0;
    for (var i = 0; i < inCore.length; i += 1) {
      var du = inCore[i][ax.u] - cx;
      var dv = inCore[i][ax.v] - cy;
      r = Math.max(r, Math.sqrt(du * du + dv * dv));
    }
    return { cx: cx, cy: cy, radius: r || 20 };
  }

  /* One leaflet: beads on a jittered hex lattice filling the annulus inner..outer. */
  function leaflet(viewer, spec, plane, fp, color) {
    if (typeof viewer.addSphere !== "function") return 0;
    var ax = AXES[opt(spec, "axis")] || AXES.z;
    var clearance = spec.clearance === undefined ? DEFAULTS.clearance : spec.clearance;
    var inner = clearance === null ? 0 : Math.max(fp.radius + clearance, 0);
    var outer = fp.radius + opt(spec, "extent");
    var target = opt(spec, "beads");
    var beadR = opt(spec, "radius");
    var alpha = opt(spec, "opacity");
    var jit = opt(spec, "jitter");
    var thick = opt(spec, "thickness");

    // Step solved from the annulus area, so the bead count stays near `target` whatever the
    // protein's girth — a wide antigen thins the lattice rather than the frame rate.
    var step = Math.max(6.8, Math.sqrt(Math.PI * (outer * outer - inner * inner) / target));
    var rowH = step * 0.866;                       // hex packing, like lipid head groups
    var jMax = Math.ceil(outer / rowH);
    var iMax = Math.ceil(outer / step) + 1;
    var drawn = 0;

    for (var j = -jMax; j <= jMax; j += 1) {
      var v0 = j * rowH;
      var offset = j % 2 === 0 ? 0 : step / 2;
      for (var i = -iMax; i <= iMax; i += 1) {
        var u0 = i * step + offset;
        var rr = Math.sqrt(u0 * u0 + v0 * v0);
        if (rr > outer || rr < inner) continue;

        // Deterministic jitter. Math.random would make the lattice shimmer on every
        // repaint, which reads as the membrane moving when only the camera did.
        var a = Math.sin(i * 12.9898 + j * 78.233) * 43758.5453; a -= Math.floor(a);
        var b = Math.sin(i * 39.3468 + j * 11.1355) * 24634.6345; b -= Math.floor(b);
        var c = Math.sin(i * 63.7264 + j * 29.8811) * 18927.1234; c -= Math.floor(c);

        var center = {};
        center[ax.u] = fp.cx + u0 + (a - 0.5) * step * jit;
        center[ax.v] = fp.cy + v0 + (c - 0.5) * step * jit;
        center[ax.n] = plane + (b - 0.5) * thick;
        viewer.addSphere({ center: center, radius: beadR, color: color, opacity: alpha });
        drawn += 1;
      }
    }
    return drawn;
  }

  /* The greasy core between the leaflets. Half-width 0.7*outer keeps the square inside the
   * bead disc so no corners protrude past it. */
  function haze(viewer, spec, fp) {
    if (!spec.haze || typeof viewer.addBox !== "function" || !spec.core) return;
    var ax = AXES[opt(spec, "axis")] || AXES.z;
    var outer = fp.radius + opt(spec, "extent");
    var s = outer * 0.7;
    var span = spec.core[1] - spec.core[0];
    var corner = {};
    corner[ax.u] = fp.cx - s;
    corner[ax.v] = fp.cy - s;
    corner[ax.n] = spec.core[0];
    // addBox dimensions are always (w, h, d) = (x, y, z), so map by axis rather than by role.
    var dims = { w: 2 * s, h: 2 * s, d: 2 * s };
    dims[{ x: "w", y: "h", z: "d" }[ax.n]] = span;
    viewer.addBox({
      corner: corner,
      dimensions: dims,
      color: spec.haze.color || "#64748b",
      opacity: spec.haze.opacity === undefined ? 0.07 : spec.haze.opacity
    });
  }

  function draw(viewer, spec) {
    spec = spec || {};
    var fp = spec.footprint || footprint(viewer, spec);
    if (!fp) return null;
    haze(viewer, spec, fp);
    var beads = 0;
    (spec.leaflets || []).forEach(function (lf) {
      if (lf && typeof lf.at === "number") {
        beads += leaflet(viewer, spec, lf.at, fp, lf.color || "#a78bfa");
      }
    });
    return { cx: fp.cx, cy: fp.cy, radius: fp.radius, beads: beads };
  }

  root.MembraneBilayer = { draw: draw, footprint: footprint, leaflet: leaflet, AXES: AXES };
}(typeof window !== "undefined" ? window : this));
