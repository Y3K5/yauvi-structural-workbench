"""Domain-aware membrane orientation from C-alpha coordinates (numpy only).

Part of memorient. Works on any predicted or experimental structure; ships no coordinates.

Why: whole-protein hydrophobic scans mis-orient multi-domain membrane proteins such as Omp85-family
barrels with periplasmic POTRA domains (BamA, TamA and related lipoproteins): the soluble domains pull the
fit and get placed inside the bilayer. This routine first finds the membrane-embedded core, fits the
membrane only to that core, and then picks the side by topology.

Method (PPM-like in spirit, much simpler; an estimate, not a measured pose):
  1. Secondary structure from C-alpha distance rules (P-SEA-like).
  2. Beta-barrel core: strands paired by C-alpha contacts; the 2-core of the strand graph with >= 8 strands
     that forms a closed ring around a common axis (low radial spread, no big angular gap).
     The ring must have an apolar lipid-facing surface (mean KD of outward residues > 0), which separates
     membrane barrels from soluble ones (GFP, lipocalins, streptavidin).
     Helical core: apolar helices (mean KD >= 1.0) of >= 14 residues that span >= 18 A.
  3. Axis: dominant eigenvector of the summed strand (or helix) direction dyads, then a restrained scan over
     normal tilt (<= 25 deg from the axis) and depth (1 A steps) maximising a burial score: outward-facing
     core residues in the slab weighted by Kyte-Doolittle hydropathy, an aromatic-girdle term for Tyr/Trp near
     both interfaces, and a penalty for non-core residues placed inside the slab.
  4. Side (sign of the normal, +normal = outside / extracellular):
       barrels  - residues before the first barrel strand (N-terminal POTRA-like domains, lipid anchor) are
                  periplasmic; longer loops face outside; the barrel N-terminus starts on the periplasmic side.
       helices  - positive-inside rule (Lys/Arg in flanking loops).
     Votes are reported; confidence is high only if the evidence agrees.
  5. Lipoproteins without a membrane core: the lipidated Cys is placed at the inner-leaflet headgroups and the
     body on the periplasmic side (side is a default, reported as low confidence).
"""
from __future__ import annotations
import numpy as np

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9,
          F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
THREE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L', LYS='K',
             MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V', MSE='M')


def _unit(v):
    v = np.asarray(v, float); n = np.linalg.norm(v, axis=-1, keepdims=True); n[n == 0] = 1; return v / n


def sse_psea(ca, resnum):
    """'a' helix, 'b' strand, 'c' coil from C-alpha distances (continuous numbering required)."""
    n = len(ca); s = np.array(['c'] * n)
    d = lambda i, j: float(np.linalg.norm(ca[i] - ca[j]))
    for i in range(n - 4):
        if resnum[i + 4] != resnum[i] + 4: continue
        d2, d3, d4 = d(i, i + 2), d(i, i + 3), d(i, i + 4)
        if abs(d2 - 5.5) < .5 and abs(d3 - 5.3) < .5 and abs(d4 - 6.4) < .6: s[i:i + 5] = 'a'
        elif abs(d2 - 6.7) < .6 and abs(d3 - 9.9) < .9 and abs(d4 - 12.4) < 1.1:
            for k in range(i, i + 5):
                if s[k] != 'a': s[k] = 'b'
    return s


def _runs(sse, kind, minlen):
    out, i, n = [], 0, len(sse)
    while i < n:
        if sse[i] == kind:
            j = i
            while j < n and (sse[j] == kind or (j + 1 < n and sse[j + 1] == kind and sse[j] == 'c' and kind == 'b')): j += 1
            if j - i >= minlen: out.append((i, j))
            i = j
        else: i += 1
    return out


def _split_turns(ca, runs, minlen=4):
    """split strand runs where the chain reverses direction (beta hairpins merged by the 5-residue windows)"""
    out = []
    for a, b in runs:
        cuts = [a]
        for i in range(a + 3, min(b - 2, len(ca) - 3)):
            u, v = ca[i] - ca[i - 3], ca[i + 3] - ca[i]
            if u @ v < 0 and i - cuts[-1] >= minlen: cuts.append(i)
        cuts.append(b)
        out += [(x, y) for x, y in zip(cuts[:-1], cuts[1:]) if y - x >= minlen]
    return out


def _dir(ca, a, b):
    return _unit(ca[b - 1] - ca[a])


def _pseudo_cb(ca):
    v = np.zeros_like(ca)
    v[1:-1] = 2 * ca[1:-1] - ca[:-2] - ca[2:]
    v[0], v[-1] = v[1], v[-2]
    return _unit(v)


def ladders(ca, resnum):
    """DSSP-like beta ladders from C-alpha only: residue pairs (i, j) closer than 5.5 A whose +/-2 neighbours are
    also paired (parallel or antiparallel register) and whose local chain is extended. Returns (strands, pair matrix)."""
    n = len(ca); D = np.linalg.norm(ca[:, None] - ca[None], axis=-1)
    cont = np.r_[np.diff(resnum) == 1, False]
    ok2 = np.zeros(n, bool)                                    # residue i has continuous i-2..i+2 and an extended local chain
    for i in range(2, n - 2):
        if cont[i - 2] and cont[i - 1] and cont[i] and cont[i + 1] and D[i - 2, i + 2] > 10.0: ok2[i] = True
    I, J = np.where((D < 5.5) & (np.abs(np.arange(n)[:, None] - np.arange(n)[None]) >= 3))
    pair = np.zeros((n, n), bool)
    for i, j in zip(I, J):
        if not (ok2[i] and ok2[j]): continue
        par = D[i - 2, j - 2] < 6.0 and D[i + 2, j + 2] < 6.0
        anti = j + 2 < n and j - 2 >= 0 and D[i - 2, j + 2] < 6.0 and D[i + 2, j - 2] < 6.0
        if par or anti: pair[i, j] = pair[j, i] = True
    inb = pair.any(1)
    for i in range(1, n - 1):                                  # fill single-residue bulges
        if not inb[i] and inb[i - 1] and inb[i + 1] and cont[i - 1] and cont[i]: inb[i] = True
    runs, i = [], 0
    while i < n:
        if inb[i]:
            j = i
            while j + 1 < n and inb[j + 1] and cont[j]: j += 1
            if j + 1 - i >= 3: runs.append((i, j + 1))
            i = j + 1
        else: i += 1
    return _split_turns(ca, runs, 3), pair


def find_barrel(ca, sse, min_strands=8, seq=None, cb=None, resnum=None):
    """Best closed (else open) beta-barrel ring from two strand sources: C-alpha ladders (robust on predicted models)
    and P-SEA strands paired by C-alpha contacts (robust on crystal structures)."""
    seq = seq or 'G' * len(ca); cb = _pseudo_cb(ca) if cb is None else cb
    resnum = np.arange(len(ca)) if resnum is None else resnum
    ls, pair = ladders(ca, resnum)
    adj1 = [set() for _ in ls]
    for i in range(len(ls)):
        for j in range(i + 1, len(ls)):
            if pair[ls[i][0]:ls[i][1], ls[j][0]:ls[j][1]].sum() >= 2: adj1[i].add(j); adj1[j].add(i)
    ps = _split_turns(ca, _runs(sse, 'b', 4)); adj2 = [set() for _ in ps]
    for i in range(len(ps)):
        A = ca[ps[i][0]:ps[i][1]]
        for j in range(i + 1, len(ps)):
            B = ca[ps[j][0]:ps[j][1]]
            if (np.linalg.norm(A[:, None] - B[None], axis=-1) < 5.6).sum() >= 3: adj2[i].add(j); adj2[j].add(i)
    cands = [c for c in (_barrel_from(ca, ls, adj1, min_strands, seq, cb), _barrel_from(ca, ps, adj2, min_strands, seq, cb)) if c]
    if not cands: return None
    return sorted(cands, key=lambda c: (not c['open'], len(c['strands']) <= 26, len(c['strands']), -c['cv']))[-1]


def _barrel_from(ca, strands, adj, min_strands, seq, cb):
    if len(strands) < min_strands: return None
    m = len(strands)
    def components(alive):
        comps, seen = [], set()
        for i in sorted(alive):
            if i in seen: continue
            st, comp = [i], []
            seen.add(i)
            while st:
                k = st.pop(); comp.append(k)
                for j in adj[k] & alive:
                    if j not in seen: seen.add(j); st.append(j)
            comps.append(sorted(comp))
        return comps

    def assess(comp, open_ok):
        segs = [strands[i] for i in comp]
        dirs = np.array([_dir(ca, a, b) for a, b in segs])
        axis = np.linalg.eigh(sum(np.outer(d, d) for d in dirs))[1][:, -1]
        keep = [k for k, d in enumerate(dirs) if abs(d @ axis) > .45]      # drop strands not running along the axis
        if len(keep) < min_strands: return None
        segs = [segs[k] for k in keep]; dirs = dirs[keep]
        axis = np.linalg.eigh(sum(np.outer(d, d) for d in dirs))[1][:, -1]
        for _ in range(4):                                                  # ring refinement: drop strands off the cylinder
            idx = np.concatenate([np.arange(a, b) for a, b in segs]); c = ca[idx].mean(0)
            sr, sz = [], []
            for a, b in segs:
                Q = ca[a:b] - c; zq = Q @ axis; rq = np.linalg.norm(Q - np.outer(zq, axis), axis=1); sr.append(rq.mean()); sz.append(zq.mean())
            med = float(np.median(sr)); good = [k for k in range(len(segs)) if abs(sr[k] - med) < max(5.0, .35 * med) and abs(sz[k]) < 14]
            if len(good) == len(segs) or len(good) < min_strands: break
            segs = [segs[k] for k in good]; dirs = dirs[good]
            axis = np.linalg.eigh(sum(np.outer(d, d) for d in dirs))[1][:, -1]
        if len(segs) < min_strands: return None
        idx = np.concatenate([np.arange(a, b) for a, b in segs]); P = ca[idx]; c = P.mean(0)
        R = (P - c) - np.outer((P - c) @ axis, axis); r = np.linalg.norm(R, axis=1)
        e1 = _unit(np.cross(axis, [1, 0, 0] if abs(axis[0]) < .9 else [0, 1, 0])); e2 = np.cross(axis, e1)
        ang = np.sort(np.degrees(np.arctan2(R @ e2, R @ e1)) % 360)
        gap = max(np.max(np.diff(ang)) if len(ang) > 1 else 360, 360 - ang[-1] + ang[0])
        cv = float(r.std() / max(r.mean(), 1e-6))
        ok = (cv < .35 and gap < 100) if not open_ok else (cv < .5 and gap < 170)
        if not (ok and 6 < r.mean() < 30): return None
        out_face = np.einsum('ij,ij->i', _unit(R), cb[idx]) > 0
        zc = (P - c) @ axis; mid = out_face & (np.abs(zc) < 10)
        face = float(np.mean([KD.get(seq[i], 0) for i in idx[mid]])) if mid.any() else -9.0
        span = float(np.ptp(zc))
        return dict(strands=segs, idx=idx, axis=axis, center=c, radius=float(r.mean()), cv=cv, gap=float(gap), ring=True, open=open_ok,
                    lipid_face_kd=round(face, 2), span=round(span, 1))

    alive = set(range(m)); changed = True
    while changed:                                   # 2-core: drop strands with < 2 partners (closed ring first)
        changed = False
        for i in list(alive):
            if len(adj[i] & alive) < 2: alive.discard(i); changed = True
    best = None
    for comp in components(alive):
        if len(comp) >= min_strands:
            cand = assess(comp, False)
            if cand and (best is None or len(cand['strands']) > len(best['strands'])): best = cand
    if best is None:                                 # open barrels (lateral gate / seam apart): whole sheet component
        for comp in components(set(range(m))):
            if len(comp) >= min_strands:
                cand = assess(comp, True)
                if cand and (best is None or len(cand['strands']) > len(best['strands'])): best = cand
    return best


def find_tm_helices(ca, sse, seq=None):
    """helices that could cross a bilayer: >= 14 residues, >= 18 A long and apolar (mean Kyte-Doolittle >= 1.0)"""
    hel = []
    for a, b in _runs(sse, 'a', 14):
        if np.linalg.norm(ca[b - 1] - ca[a]) < 18: continue
        if seq is not None and np.mean([KD.get(x, 0) for x in seq[a:b]]) < 1.0: continue
        hel.append((a, b))
    return hel


def _cone(axis, max_deg, step_deg):
    axis = _unit(axis); e1 = _unit(np.cross(axis, [1, 0, 0] if abs(axis[0]) < .9 else [0, 1, 0])); e2 = np.cross(axis, e1)
    out = [axis]
    for t in np.arange(step_deg, max_deg + 1e-6, step_deg):
        k = max(6, int(round(2 * np.pi * np.sin(np.radians(t)) / np.radians(step_deg))))
        for p in np.linspace(0, 2 * np.pi, k, endpoint=False):
            out.append(np.cos(np.radians(t)) * axis + np.sin(np.radians(t)) * (np.cos(p) * e1 + np.sin(p) * e2))
    return np.array(out)


def _scan(ca, seq, core_idx, axis, center, half, outward_only, cb):
    kd = np.array([KD.get(a, 0.0) for a in seq]) / 4.5
    arom = np.array([a in 'WY' for a in seq], float)
    core = np.zeros(len(ca), bool); core[core_idx] = True
    best = (-1e9, None, None)
    for n in _cone(axis, 25, 3):
        z0 = (ca - center) @ n
        if outward_only:
            rad = (ca - center) - np.outer(z0, n); facing = (np.einsum('ij,ij->i', _unit(rad), cb) > 0.0)
        else:
            facing = np.ones(len(ca), bool)
        wcore = core & facing
        for d in np.arange(-14, 14.01, 1.0):
            z = z0 - d; inside = 1 / (1 + np.exp((np.abs(z) - half) / 1.5))
            girdle = np.exp(-((np.abs(z) - (half - 1.5)) / 3.0) ** 2)
            s = (kd * inside)[wcore].sum() + 1.2 * (arom * girdle)[core].sum() - 0.6 * inside[~core].sum()
            if s > best[0]: best = (s, n, d)
    s, n, d = best
    return n, center + d * n, float(s)


def _frame(n):
    e1 = _unit(np.cross(n, [1, 0, 0] if abs(n[0]) < .9 else [0, 1, 0])); e2 = np.cross(e1, n)
    return np.array([e1, n, e2])           # rows: x, y (= membrane normal, +y outside), z


def orient(ca, seq, resnum=None, lipobox_index=None, half=None, context='outer_membrane', skip_n=0):
    """ca: (N,3) C-alpha coordinates; seq: one-letter codes aligned to ca; resnum: residue numbers.
    lipobox_index: index (into ca) of a lipidated Cys, if known. skip_n: leading residues to ignore for the side vote
    (e.g. a predicted signal peptide, which is cleaved and floats freely in predicted models). Returns a dict."""
    ca = np.asarray(ca, float); seq = ''.join(seq); n = len(ca)
    resnum = np.arange(1, n + 1) if resnum is None else np.asarray(resnum)
    sse = sse_psea(ca, resnum); cb = _pseudo_cb(ca)
    out = dict(method=None, n_core=0, votes=[], confidence='low')
    bar = find_barrel(ca, sse, seq=seq, cb=cb, resnum=resnum)
    if bar:                                                   # lipid-facing hydropathy at the best slab placement
        h0 = 12.0 if half is None else half
        n0, o0, _ = _scan(ca, seq, bar['idx'], bar['axis'], bar['center'], h0, True, cb)
        P = ca[bar['idx']]; z0 = (P - o0) @ n0; Rr = (P - o0) - np.outer(z0, n0)
        face = np.einsum('ij,ij->i', _unit(Rr), cb[bar['idx']]) > 0; m = face & (np.abs(z0) < h0 - 2)
        bar['lipid_face_kd'] = round(float(np.mean([KD.get(seq[i], 0) for i in bar['idx'][m]])) if m.any() else -9.0, 2)
        if bar['lipid_face_kd'] <= (0.3 if bar['open'] else 0.0):  # polar outer face: a soluble barrel, not a membrane one
            out['rejected_barrel'] = dict(strands=len(bar['strands']), lipid_face_kd=bar['lipid_face_kd'], open=bar['open']); bar = None
    if bar:
        h = 12.0 if half is None else half
        nrm, origin, score = _scan(ca, seq, bar['idx'], bar['axis'], bar['center'], h, True, cb)
        first = int(min(bar['idx'])); last = int(max(bar['idx'])) + 1
        z = (ca - origin) @ nrm; votes = []                  # (rule, side that is OUTSIDE (+1 = +normal), value, weight)
        pre = np.arange(min(skip_n, first), first); pz = z[pre]; up, dn = int((pz > h + 8).sum()), int((pz < -h - 8).sum())   # well clear of the slab: plug loops near the slab do not count
        if up + dn >= 20 and max(up, dn) >= .65 * (up + dn):  # N-terminal domains outside the slab (POTRA-like, TonB box) -> periplasm
            votes.append(('n_terminal_domains_periplasmic', -1 if up > dn else 1, [up, dn], 2))   # plugs inside the barrel lumen sit in the slab and do not vote
        loops = np.setdiff1d(np.arange(first, last), bar['idx'])
        up, dn = int((z[loops] > h).sum()), int((z[loops] < -h).sum())
        if abs(up - dn) >= 6: votes.append(('longer_loops_outside', 1 if up > dn else -1, [up, dn], 2 if max(up, dn) >= 3 * max(1, min(up, dn)) and abs(up - dn) >= 20 else 1))
        zs = float(z[bar['strands'][0][0]])
        if abs(zs) > 4: votes.append(('first_strand_starts_periplasmic', 1 if zs < 0 else -1, round(zs, 1), 1))
        if lipobox_index is not None and lipobox_index < len(z) and abs(z[lipobox_index]) > 4:
            votes.append(('lipid_anchor_periplasmic_leaflet', 1 if z[lipobox_index] < 0 else -1, round(float(z[lipobox_index]), 1), 1))
        tot = sum(v[1] * v[3] for v in votes); sign = 1 if tot >= 0 else -1
        agree = all(v[1] == sign for v in votes)
        nrm = nrm * sign
        tilt = float(np.degrees(np.arccos(min(1, abs(nrm @ bar['axis'])))))
        out.update(method='beta_barrel_core' + ('_open' if bar['open'] else ''), n_core=len(bar['strands']), radius=round(bar['radius'], 1), ring_cv=round(bar['cv'], 2),
                   votes=[dict(rule=v[0], says='as_chosen' if v[1] == sign else 'against', value=v[2], weight=v[3]) for v in votes],
                   confidence='high' if agree and len(votes) >= 2 and (bar['cv'] < .25 if bar['open'] else bar['cv'] < .32) else 'medium' if agree or abs(tot) >= 3 else 'low',
                   score=round(score, 2), lipid_face_kd=bar['lipid_face_kd'], half=h, tilt_vs_core_axis=round(tilt, 1),
                   depth_offset=round(float((origin - bar['center']) @ nrm), 1), core_idx=bar['idx'].tolist())
    else:
        hel = [hh for hh in find_tm_helices(ca, sse, seq) if hh[0] >= skip_n]   # a predicted signal peptide is not a TM helix
        if hel and context != 'soluble':
            h = 15.0 if half is None else half
            dirs = np.array([_dir(ca, a, b) for a, b in hel]); M = sum(np.outer(d, d) for d in dirs); axis = np.linalg.eigh(M)[1][:, -1]
            idx = np.concatenate([np.arange(a, b) for a, b in hel]); c = ca[idx].mean(0)
            nrm, origin, score = _scan(ca, seq, idx, axis, c, h, False, cb)
            z = (ca - origin) @ nrm; core = np.zeros(n, bool); core[idx] = True
            flank = (~core) & (np.abs(z) > h - 2) & (np.abs(z) < h + 12)
            pos = np.array([a in 'KR' for a in seq])
            up, dn = int((pos & flank & (z > 0)).sum()), int((pos & flank & (z < 0)).sum())
            votes = []
            if up != dn: votes.append(('positive_inside', -1 if up > dn else 1, [up, dn]))
            sign = votes[0][1] if votes else 1; nrm = nrm * sign
            single = len(hel) == 1
            out.update(method='helical_core', n_core=len(hel), half=h, score=round(score, 2),
                       votes=[dict(rule=v[0], says='as_chosen', value=v[2]) for v in votes],
                       confidence='low' if single or not votes or abs(up - dn) < 2 else 'medium',
                       tilt_vs_core_axis=round(float(np.degrees(np.arccos(min(1, abs(nrm @ axis))))), 1),
                       depth_offset=round(float((origin - c) @ nrm), 1), core_idx=idx.tolist())
        elif lipobox_index is not None and lipobox_index < n:
            h = 12.0 if half is None else half
            body = ca.mean(0); nrm = _unit(ca[lipobox_index] - body)      # anchor above the body
            origin = ca[lipobox_index] + (h - 2.0) * nrm                   # Cys sits in the inner-leaflet headgroups
            out.update(method='lipid_anchor', n_core=0, half=h, votes=[dict(rule='default_periplasmic_body', says='as_chosen', value=None)],
                       confidence='low', tilt_vs_core_axis=None, depth_offset=None, core_idx=[])
        else:
            out.update(method='none'); return out
    R = _frame(nrm)
    out.update(normal=[round(float(x), 4) for x in nrm], origin=[round(float(x), 2) for x in origin],
               rot=[round(float(x), 5) for x in R.ravel()])
    zz = (ca - origin) @ nrm; out['in_slab'] = int((np.abs(zz) <= out['half']).sum())
    return out


def read_pdb_ca(path, chain=None):
    ca, seq, num = [], [], []; seen = set()
    for ln in open(path):
        if ln.startswith('ENDMDL'): break
        if not ln.startswith(('ATOM', 'HETATM')) or ln[12:16].strip() != 'CA': continue
        if chain and ln[21] != chain: continue
        rn = ln[17:20]; key = (ln[21], ln[22:27])
        if rn not in THREE or key in seen: continue
        if ln[16] not in ' A': continue
        seen.add(key); ca.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]); seq.append(THREE[rn]); num.append(int(ln[22:26]))
    return np.array(ca), ''.join(seq), np.array(num)


def main(argv=None):
    """python -m memorient.domain_orient STRUCTURE.pdb [--chain A] [--half 12] [--json out.json]"""
    import argparse, json
    ap = argparse.ArgumentParser(description='Domain-aware membrane orientation (C-alpha only, offline).')
    ap.add_argument('pdb'); ap.add_argument('--chain'); ap.add_argument('--half', type=float); ap.add_argument('--skip-n', type=int, default=0)
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    ca, seq, num = read_pdb_ca(a.pdb, a.chain)
    res = orient(ca, seq, num, half=a.half, skip_n=a.skip_n); res.pop('core_idx', None)
    txt = json.dumps(res, indent=1)
    if a.json: open(a.json, 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
