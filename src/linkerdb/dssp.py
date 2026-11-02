"""Minimal DSSP (Kabsch & Sander, 1983) secondary-structure assignment in numpy.

Implements the parts of DSSP needed to decide whether a residue is in *regular*
secondary structure:

* electrostatic H-bond energy ``E = 0.084*332*(1/rON + 1/rCH - 1/rOH - 1/rCN)``,
  bond if ``E < -0.5 kcal/mol``; amide H placed from the previous C=O;
  like DSSP, each NH keeps only its two best acceptors and Pro is no donor;
* n-turns (n = 3, 4, 5) and minimal helices -> ``H`` (alpha), ``G`` (3-10), ``I`` (pi);
* parallel / antiparallel bridges -> ``E`` (ladder of >= 2 consecutive bridges) or ``B``.

Not implemented: beta-bulge merging of ladders, bends/turn letters (``T``, ``S``).
Unobserved positions are ``'?'``. Priority: H > E > B > G > I (DSSP 2.x).
"""

from __future__ import annotations

import numpy as np

HBOND_E = -0.5
Q = 0.084 * 332.0
CA_CUTOFF = 9.0
PEPTIDE_MAX = 2.5  # DSSP's chain-break C-N distance
REGULAR = frozenset("HEI")


def _amide_h(n: np.ndarray, c: np.ndarray, o: np.ndarray, linked_prev: np.ndarray) -> np.ndarray:
    h = np.full_like(n, np.nan)
    co = c[:-1] - o[:-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        co /= np.linalg.norm(co, axis=1, keepdims=True)
    h[1:] = n[1:] + co
    h[~linked_prev] = np.nan
    return h


def hbond_matrix(coords: np.ndarray, is_pro: np.ndarray | None = None, block: int = 512) -> np.ndarray:
    """Boolean ``hb[d, a]``: NH of residue ``d`` donates to C=O of residue ``a``."""
    n_res = len(coords)
    n, ca, c, o = (coords[:, k] for k in range(4))
    link = links(coords)
    linked_prev = np.zeros(n_res, dtype=bool)
    linked_prev[1:] = link
    h = _amide_h(n, c, o, linked_prev)
    donor_ok = ~np.isnan(h).any(axis=1)
    if is_pro is not None:
        donor_ok &= ~np.asarray(is_pro, dtype=bool)
    acc_ok = ~np.isnan(c).any(axis=1) & ~np.isnan(o).any(axis=1)
    hb = np.zeros((n_res, n_res), dtype=bool)
    idx = np.arange(n_res)
    for start in range(0, n_res, block):
        d = idx[start:start + block]
        with np.errstate(invalid="ignore", divide="ignore"):
            dca = np.linalg.norm(ca[d, None, :] - ca[None, :, :], axis=2)
            r_on = np.linalg.norm(o[None, :, :] - n[d, None, :], axis=2)
            r_ch = np.linalg.norm(c[None, :, :] - h[d, None, :], axis=2)
            r_oh = np.linalg.norm(o[None, :, :] - h[d, None, :], axis=2)
            r_cn = np.linalg.norm(c[None, :, :] - n[d, None, :], axis=2)
            e = Q * (1 / r_on + 1 / r_ch - 1 / r_oh - 1 / r_cn)
        cand = (dca < CA_CUTOFF) & donor_ok[d, None] & acc_ok[None, :]
        cand &= np.abs(d[:, None] - idx[None, :]) >= 2
        e = np.where(cand & ~np.isnan(e), e, np.inf)
        e = np.maximum(e, -9.9)  # DSSP clamps the minimum energy
        # keep the two best acceptors per donor (as DSSP does)
        if n_res > 2:
            best = np.argpartition(e, 1, axis=1)[:, :2]
        else:
            best = np.argsort(e, axis=1)[:, :2]
        rows = np.repeat(np.arange(len(d)), best.shape[1])
        cols = best.ravel()
        ok = e[rows, cols] < HBOND_E
        hb[d[rows[ok]], cols[ok]] = True
    return hb


def links(coords: np.ndarray) -> np.ndarray:
    """``links[i]`` True if residues ``i`` and ``i+1`` are peptide-bonded (no chain break)."""
    if len(coords) < 2:
        return np.zeros(0, dtype=bool)
    with np.errstate(invalid="ignore"):
        d = np.linalg.norm(coords[:-1, 2] - coords[1:, 0], axis=1)
        return d < PEPTIDE_MAX


def assign(coords: np.ndarray, is_pro: np.ndarray | None = None) -> str:
    """DSSP-like secondary structure string for an (N, 4, 3) N/CA/C/O array."""
    n_res = len(coords)
    observed = ~np.isnan(coords[:, 1]).any(axis=1)
    ss = np.array(["-"] * n_res, dtype="<U1")
    if n_res < 3:
        ss[~observed] = "?"
        return "".join(ss)
    hb = hbond_matrix(coords, is_pro)
    link = links(coords)
    # no_break[i, j] for j > i: all links in i..j-1 present (cumulative count trick)
    cum = np.concatenate([[0], np.cumsum(~link)])

    def no_break(i: int, j: int) -> bool:
        return 0 <= i and j < n_res and cum[j] - cum[i] == 0

    def bond(d: int, a: int) -> bool:
        return 0 <= d < n_res and 0 <= a < n_res and hb[d, a]

    # turns and helices ---------------------------------------------------
    helix = {}
    for stride, letter in ((4, "H"), (3, "G"), (5, "I")):
        turn = np.zeros(n_res, dtype=bool)
        for i in range(n_res - stride):
            if no_break(i, i + stride) and hb[i + stride, i]:
                turn[i] = True
        mark = np.zeros(n_res, dtype=bool)
        for i in range(1, n_res - stride):
            if turn[i - 1] and turn[i]:
                mark[i:i + stride] = True
        helix[letter] = mark

    # bridges -------------------------------------------------------------
    par = []
    anti = []
    cand_i, cand_j = np.nonzero(np.triu(hb | hb.T, 1))
    pairs = set()
    for a, b in zip(cand_i, cand_j):
        # any residue pair (i, j) whose bridge test involves bond a<->b lies within +-1
        for i in (a - 1, a, a + 1):
            for j in (b - 1, b, b + 1):
                if 1 <= i < j - 2 and j < n_res - 1:
                    pairs.add((min(i, j), max(i, j)))
    for i, j in sorted(pairs):
        if not (no_break(i - 1, i + 1) and no_break(j - 1, j + 1)):
            continue
        if (bond(i + 1, j) and bond(j, i - 1)) or (bond(j + 1, i) and bond(i, j - 1)):
            par.append((i, j))
        elif (bond(i + 1, j - 1) and bond(j + 1, i - 1)) or (bond(j, i) and bond(i, j)):
            anti.append((i, j))
    strand, bridge = _ladders(par, anti, n_res, no_break)

    # priority: lowest first, overwritten by higher ------------------------
    ss[helix["I"]] = "I"
    ss[helix["G"]] = "G"
    ss[bridge] = "B"
    ss[strand] = "E"
    ss[helix["H"]] = "H"
    ss[~observed] = "?"
    return "".join(ss)


def _ladders(par, anti, n_res, no_break):
    """Group bridges into ladders, merge bulge-linked ladders (DSSP rules) and
    return boolean arrays ``(strand, isolated_bridge)``."""
    bridges = sorted([(i, j, "p") for i, j in par] + [(i, j, "a") for i, j in anti])
    ladders: list[dict] = []
    for i, j, t in bridges:
        for lad in ladders:
            li, lj = lad["pairs"][-1]
            if lad["type"] != t or i != li + 1 or not no_break(li, i):
                continue
            if (t == "p" and j == lj + 1 and no_break(lj, j)) or (t == "a" and j == lj - 1 and no_break(j, lj)):
                lad["pairs"].append((i, j))
                break
        else:
            ladders.append({"type": t, "pairs": [(i, j)], "linked": False})
    # bulge merging: union-find over ladders
    parent = list(range(len(ladders)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    span = []
    for lad in ladders:
        iis = [p[0] for p in lad["pairs"]]
        jjs = [p[1] for p in lad["pairs"]]
        span.append((min(iis), max(iis), min(jjs), max(jjs)))
    for a in range(len(ladders)):
        ibi, iei, jbi, jei = span[a]
        for b in range(a + 1, len(ladders)):
            ibj, iej, jbj, jej = span[b]
            if ladders[a]["type"] != ladders[b]["type"]:
                continue
            if not no_break(min(ibi, ibj), max(iei, iej)) or not no_break(min(jbi, jbj), max(jei, jej)):
                continue
            # DSSP computes these gaps in unsigned arithmetic: negative gaps never qualify
            d_i = ibj - iei
            if not 0 <= d_i < 6 or (iei >= ibj and ibi <= iej):
                continue
            d_j = (jbj - jei) if ladders[a]["type"] == "p" else (jbi - jej)
            bulge = (0 <= d_j < 6 and d_i < 3) or (0 <= d_j < 3)
            if bulge:
                ladders[a]["linked"] = ladders[b]["linked"] = True
                parent[find(b)] = find(a)
    strand = np.zeros(n_res, dtype=bool)
    bridge = np.zeros(n_res, dtype=bool)
    groups: dict[int, list[int]] = {}
    for k in range(len(ladders)):
        groups.setdefault(find(k), []).append(k)
    for members in groups.values():
        n_pairs = sum(len(ladders[k]["pairs"]) for k in members)
        linked = len(members) > 1
        for k in members:
            ibi, iei, jbi, jei = span[k]
            if n_pairs > 1 or linked:
                strand[ibi:iei + 1] = True
                strand[jbi:jei + 1] = True
            else:
                bridge[ibi] = bridge[jbi] = True
        if linked:  # fill the bulge residues between linked ladders
            iis = [x for k in members for x in span[k][:2]]
            jjs = [x for k in members for x in span[k][2:]]
            strand[min(iis):max(iis) + 1] = True
            strand[min(jjs):max(jjs) + 1] = True
    return strand, bridge & ~strand


def is_regular(code: str) -> bool:
    return code in REGULAR
