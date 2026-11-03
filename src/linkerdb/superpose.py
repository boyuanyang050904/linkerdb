"""Superposition-based inter-domain mobility — a rotation-aware rigidity measure.

The dom_dist / elbow-angle rigidity call (aggregate.py) is a cheap proxy. The clean
measure of how two domains move relative to each other is: superpose the N-domain
across a protein's structures, then see how far the C-domain has moved. This is
computed here.

For two structures of the same linker:

* correspond residues by **UniProt position** (SEQRES numbering differs per entry);
* Kabsch-superpose the N-domain Cα atoms (shared positions);
* apply that transform to the C-domain and measure its Cα RMSD vs the reference
  (``c_rmsd`` — total inter-domain displacement, captures translation *and* rotation);
* report the inter-domain rotation angle between the two domains' frames (``hinge``).

A group's mobility is the max over its structures against the best-resolution
reference. It still cannot see flexibility the PDB does not sample (if every
structure is the same conformer, mobility is ~0 — a data limit, not a metric one).

This is heavier than the proxy (it re-reads structures and re-derives domains), so it
is an analysis/validation tool here, not part of the default build.
"""

from __future__ import annotations

import logging

import gemmi
import numpy as np


log = logging.getLogger(__name__)

RIGID_RMSD = 1.5     # Å: C-domain displacement below this (after N-superposition) = rigid
FLEX_RMSD = 3.0      # Å: above this = flexible
RIGID_HINGE = 10.0   # degrees
FLEX_HINGE = 20.0
MIN_CORE = 5         # min shared Cα to superpose a domain


def domain_unp_ca(chain, mapped_domain, uniprot: str) -> dict[int, np.ndarray]:
    """``{uniprot_position: CA xyz}`` for modelled residues of this domain that map to ``uniprot``."""
    out: dict[int, np.ndarray] = {}
    for s, e in mapped_domain.segments:
        for sid in range(s, e + 1):
            if not (1 <= sid <= chain.n) or not chain.res(sid).observed:
                continue
            hit = chain.unp_position(sid)
            if hit is None or hit[0] != uniprot:
                continue
            p = chain.coords[sid - 1, 1]
            if not np.isnan(p).any():
                out[hit[1]] = p
    return out


def _rot_angle(R: np.ndarray) -> float:
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1.0, 1.0))))


def _positions(d: dict[int, np.ndarray], keys) -> list:
    return [gemmi.Position(*d[k]) for k in keys]


def relative_shift(ref_n, ref_c, mov_n, mov_c) -> tuple[float, float] | None:
    """C-domain RMSD and inter-domain hinge angle of ``mov`` relative to ``ref``."""
    n_keys = sorted(set(ref_n) & set(mov_n))
    c_keys = sorted(set(ref_c) & set(mov_c))
    if len(n_keys) < MIN_CORE or len(c_keys) < MIN_CORE:
        return None
    sup_n = gemmi.superpose_positions(_positions(ref_n, n_keys), _positions(mov_n, n_keys))
    Rn, tn = np.array(sup_n.transform.mat.tolist()), np.array(sup_n.transform.vec.tolist())
    ref_c_arr = np.array([ref_c[k] for k in c_keys])
    mov_c_arr = np.array([mov_c[k] for k in c_keys])
    moved = mov_c_arr @ Rn.T + tn                       # C of mov, in ref's N-frame
    c_rmsd = float(np.sqrt(((moved - ref_c_arr) ** 2).sum(axis=1).mean()))
    sup_c = gemmi.superpose_positions(_positions(ref_c, c_keys), _positions(mov_c, c_keys))
    Rc = np.array(sup_c.transform.mat.tolist())
    hinge = _rot_angle(Rn.T @ Rc)   # rotation of the C-domain frame relative to the N-domain frame
    return c_rmsd, hinge


def mobility(domains_by_structure: list[tuple[str, dict, dict]]) -> dict:
    """``domains_by_structure``: list of (label, n_domain_ca, c_domain_ca), reference first.

    Returns max C-domain RMSD and max hinge angle of the rest against the reference.
    """
    if len(domains_by_structure) < 2:
        return {"n_used": len(domains_by_structure), "max_c_rmsd": None, "max_hinge": None, "n_pairs": 0}
    _, ref_n, ref_c = domains_by_structure[0]
    rmsds, hinges, n = [], [], 0
    for _, mn, mc in domains_by_structure[1:]:
        res = relative_shift(ref_n, ref_c, mn, mc)
        if res is None:
            continue
        rmsds.append(res[0])
        hinges.append(res[1])
        n += 1
    if not n:
        return {"n_used": len(domains_by_structure), "max_c_rmsd": None, "max_hinge": None, "n_pairs": 0}
    return {"n_used": len(domains_by_structure), "max_c_rmsd": max(rmsds), "max_hinge": max(hinges), "n_pairs": n}


def classify(max_c_rmsd: float | None, max_hinge: float | None) -> str:
    if max_c_rmsd is None:
        return "undetermined"
    if max_c_rmsd <= RIGID_RMSD and (max_hinge is None or max_hinge <= RIGID_HINGE):
        return "rigid"
    if max_c_rmsd >= FLEX_RMSD or (max_hinge is not None and max_hinge >= FLEX_HINGE):
        return "flexible"
    return "intermediate"
