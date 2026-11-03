"""Second layer: group linker *instances* into *canonical linkers*.

An instance is one domain junction seen in one chain. The same natural linker
appears in many PDB structures; grouping those instances is what lets us say
whether a linker is rigid or flexible (how much its geometry varies across
structures) — information the per-structure non-redundant set throws away.

Grouping key: UniProt accession + overlapping UniProt interval. Instances of the
same junction across structures share a UniProt region and merge; distinct
junctions sit at different regions and stay separate. Instances without a
UniProt mapping cannot be grouped across structures and are emitted as singletons
flagged ``ungrouped``.

Rigidity is judged from the spread of the end-to-end distance and elbow angle
across structures (needs >= ``min_struct_for_rigidity`` structures); otherwise it
falls back to the consensus secondary structure with low confidence. The
thresholds here are heuristic and are the first thing to recalibrate against a
hand-curated golden set (see docs/design.md roadmap).
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from . import developability

# rigidity thresholds (coefficient of variation of domain-centroid distance; std of elbow angle in degrees)
DD_CV_RIGID = 0.05
DD_CV_FLEX = 0.15
ANGLE_STD_RIGID = 15.0
ANGLE_STD_FLEX = 40.0
MIN_STRUCT_FOR_RIGIDITY = 3


def _merge_intervals(rows: list[tuple[int, int, int]]) -> list[list[int]]:
    """Union-find over (lo, hi, idx); return groups of idx with overlapping [lo, hi]."""
    rows = sorted(rows)
    groups: list[list[int]] = []
    cur: list[int] = []
    cur_hi = None
    for lo, hi, idx in rows:
        if cur and lo <= cur_hi:
            cur.append(idx)
            cur_hi = max(cur_hi, hi)
        else:
            if cur:
                groups.append(cur)
            cur, cur_hi = [idx], hi
    if cur:
        groups.append(cur)
    return groups


def _modal(series):
    vals = [v for v in series if isinstance(v, str) and v]
    return Counter(vals).most_common(1)[0][0] if vals else None


def _rigidity(g: pd.DataFrame) -> tuple[str, str, float | None, float | None]:
    """Return (label, confidence, dom_dist_cv, angle_std).

    Rigidity is judged from the spread of the **domain-centroid distance** across
    structures, not the end-to-end distance: the latter depends on where the linker
    boundary is drawn, which CATH/DSSP place differently in different structures, so
    it carries boundary noise. dom_dist uses the whole flanking domains and is
    boundary-independent — it directly measures whether the two domains stay a fixed
    distance apart. (It cannot see a pure hinge rotation at constant COM distance; see
    design.md known limitations.)
    """
    dd = g["dom_dist"].dropna().to_numpy(dtype=float)
    ang = g["elbow_angle"].dropna().to_numpy(dtype=float)
    frac_helix = float((g["junction_ss"] == "helix_spanning").mean())
    frac_disorder = float((g["n_missing"] > 0).mean())
    dd_cv = float(dd.std() / dd.mean()) if len(dd) >= 2 and dd.mean() > 0 else None
    angle_std = float(ang.std()) if len(ang) >= 2 else None
    n = g[["pdb_id", "auth_asym_id"]].drop_duplicates().shape[0]

    if frac_helix >= 0.5:
        return "rigid_helical", ("high" if n >= MIN_STRUCT_FOR_RIGIDITY else "low"), dd_cv, angle_std
    if n >= MIN_STRUCT_FOR_RIGIDITY and len(dd) >= MIN_STRUCT_FOR_RIGIDITY:
        rigid = (dd_cv is not None and dd_cv <= DD_CV_RIGID) and (angle_std is None or angle_std <= ANGLE_STD_RIGID)
        flex = (dd_cv is not None and dd_cv >= DD_CV_FLEX) or (angle_std is not None and angle_std >= ANGLE_STD_FLEX)
        label = "rigid" if rigid else ("flexible" if flex else "intermediate")
        return label, "high", dd_cv, angle_std
    if frac_disorder >= 0.5:
        return "flexible_disordered", "low", dd_cv, angle_std
    # single/few structures: fall back to the dominant secondary structure
    ss = _modal(g["junction_ss"])
    guess = {"helix_spanning": "rigid_helical", "strand_spanning": "rigid",
             "coil": "flexible", "sse_abutting": "intermediate"}.get(ss, "undetermined")
    return guess, "low", dd_cv, angle_std


def _canonical_row(cid: str, g: pd.DataFrame, grouped: bool) -> dict:
    # representative instance: most observed, then best resolution
    g = g.assign(_obs=g["ss_length"] - g["n_missing"])
    rep = g.sort_values(["_obs", "resolution"], ascending=[False, True]).iloc[0]
    n_struct = g[["pdb_id", "auth_asym_id"]].drop_duplicates().shape[0]
    label, conf, dd_cv, angle_std = _rigidity(g)
    seq = rep["ss_sequence"] or ""
    liab = developability.scan(seq)
    e2e = g["end_to_end"].dropna().to_numpy(dtype=float)
    dd = g["dom_dist"].dropna().to_numpy(dtype=float)
    ang = g["elbow_angle"].dropna().to_numpy(dtype=float)
    clusters = sorted({c for c in g.get("cluster40", pd.Series(dtype=object)).dropna()})
    return {
        "canonical_id": cid,
        "uniprot": rep["uniprot"],
        "unp_lo": int(g["unp_lo"].median()) if grouped else None,
        "unp_hi": int(g["unp_hi"].median()) if grouped else None,
        "grouped": grouped,
        "n_structures": n_struct,
        "n_entries": g["pdb_id"].nunique(),
        "rep_linker_id": rep["linker_id"],
        "rep_pdb_id": rep["pdb_id"],
        "rep_chain": rep["auth_asym_id"],
        "best_resolution": float(g["resolution"].min()) if g["resolution"].notna().any() else None,
        "sequence": seq,
        "length_median": float(g["ss_length"].median()),
        "length_min": int(g["ss_length"].min()),
        "length_max": int(g["ss_length"].max()),
        "rigidity": label,
        "rigidity_confidence": conf,
        "end_to_end_median": float(np.median(e2e)) if len(e2e) else None,
        "end_to_end_cv": float(e2e.std() / e2e.mean()) if len(e2e) >= 2 and e2e.mean() > 0 else None,
        "elbow_angle_median": float(np.median(ang)) if len(ang) else None,
        "elbow_angle_std": angle_std,
        "dom_dist_median": float(np.median(dd)) if len(dd) else None,
        "dom_dist_cv": dd_cv,
        "frac_disordered": float((g["n_missing"] > 0).mean()),
        "frac_helix": float((g["junction_ss"] == "helix_spanning").mean()),
        "junction_ss": _modal(g["junction_ss"]),
        "junction_type": _modal(g["junction_type"]),
        "cath_n": _modal(g["cath_n"]),
        "cath_c": _modal(g["cath_c"]),
        "taxid": _modal(g["taxid"].astype(str)),
        "is_human": bool((g["taxid"].astype(str) == "9606").any()),
        "description": _modal(g["description"]) if "description" in g else None,
        "n_clusters40": len(clusters),
        **liab,
        "liability_tags": developability.summary_flags(liab),
        "liability_score": developability.liability_score(liab),
    }


def aggregate(instances: pd.DataFrame, qc_only: bool = True) -> pd.DataFrame:
    """Group linker instances into canonical linkers (one row each)."""
    df = instances[instances["qc_pass"]] if qc_only and "qc_pass" in instances else instances
    df = df[df["ss_length"] > 0] if "ss_length" in df else df   # a linker with no residues is not a sample
    rows: list[dict] = []
    # grouped: have a UniProt accession and interval
    mapped = df[df["uniprot"].notna() & df["unp_lo"].notna() & df["unp_hi"].notna()]
    for acc, g in mapped.groupby("uniprot"):
        g = g.reset_index(drop=True)
        intervals = [(int(min(lo, hi)), int(max(lo, hi)), i)
                     for i, (lo, hi) in enumerate(zip(g["unp_lo"], g["unp_hi"]))]
        for members in _merge_intervals(intervals):
            sub = g.iloc[members]
            lo, hi = int(sub["unp_lo"].median()), int(sub["unp_hi"].median())
            rows.append(_canonical_row(f"{acc}:{lo}-{hi}", sub, grouped=True))
    # ungrouped: no usable UniProt mapping -> one canonical per instance
    for _, r in df[~(df["uniprot"].notna() & df["unp_lo"].notna() & df["unp_hi"].notna())].iterrows():
        rows.append(_canonical_row(r["linker_id"], pd.DataFrame([r]), grouped=False)
                    | {"grouped": False})
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["n_structures", "uniprot", "unp_lo"], ascending=[False, True, True]).reset_index(
            drop=True)
    return out
