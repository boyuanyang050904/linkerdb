"""Query the canonical-linker table the way a fusion-protein designer would:
'give me natural linkers that hold two domains ~N Å apart, rigid, human, with no
glycosylation site'.

This narrows and ranks candidates; it does not predict that a linker will work in
a particular construct. Treat the output as a shortlist to model and test.
"""

from __future__ import annotations

import pandas as pd

RIGIDITY_GROUPS = {
    "rigid": {"rigid", "rigid_helical", "strand_spanning"},
    "flexible": {"flexible", "flexible_disordered"},
    "any": None,
}


def query_canonical(df: pd.DataFrame, distance: float | None = None, distance_tol: float = 5.0,
                    rigidity: str | None = None, human_only: bool = False,
                    min_length: int | None = None, max_length: int | None = None,
                    min_structures: int = 1, no_glyc: bool = False, no_cys: bool = False,
                    max_liability: int | None = None, linker_class: str | None = None,
                    top: int = 20) -> pd.DataFrame:
    """Filter and rank canonical linkers. ``distance`` is the target domain-COM distance (Å).

    Rigidity is read from ``rigidity_final`` (superposition where available, else the
    proxy) when that column is present, otherwise from the proxy ``rigidity``.
    """
    q = df.copy()
    rig_col = "rigidity_final" if "rigidity_final" in q.columns else "rigidity"
    if rigidity and rigidity != "any":
        allowed = RIGIDITY_GROUPS.get(rigidity, {rigidity})
        q = q[q[rig_col].isin(allowed)]
    if linker_class and "linker_class" in q:
        q = q[q["linker_class"] == linker_class]
    if human_only:
        q = q[q["is_human"]]
    if min_length is not None:
        q = q[q["length_median"] >= min_length]
    if max_length is not None:
        q = q[q["length_median"] <= max_length]
    if min_structures > 1:
        q = q[q["n_structures"] >= min_structures]
    if no_glyc:
        q = q[q["n_glyc"] == 0]
    if no_cys:
        q = q[q["n_cys"] == 0]
    if max_liability is not None:
        q = q[q["liability_score"] <= max_liability]
    if distance is not None:
        q = q[q["dom_dist_median"].notna()]
        q = q[(q["dom_dist_median"] - distance).abs() <= distance_tol]
        q = q.assign(_dist_err=(q["dom_dist_median"] - distance).abs())
        q = q.sort_values(["_dist_err", "n_structures", "liability_score", "best_resolution"],
                          ascending=[True, False, True, True]).drop(columns="_dist_err")
    else:
        q = q.sort_values(["n_structures", "liability_score", "best_resolution"],
                          ascending=[False, True, True])
    return q.head(top).reset_index(drop=True)


def format_hits(hits: pd.DataFrame) -> str:
    """One line per hit. ``domΔcv`` is the across-structure spread of the domain
    distance that drives the rigidity call (low = rigid); it is ``-`` when the
    linker has too few structures to measure, in which case trust ``rigidity`` less."""
    if hits.empty:
        return "no linkers match"
    lines = []
    for r in hits.itertuples():
        raw = getattr(r, "liability_tags", None)
        tags = ",".join(list(raw)) if raw is not None and len(raw) else "-"
        cls = getattr(r, "linker_class", "coil")
        rig = getattr(r, "rigidity_final", None) or r.rigidity
        src = (getattr(r, "rigidity_source", "") or "")[:4]
        lines.append(
            f"{r.canonical_id:22s} len={r.length_median:>4.0f} {cls:14s} {rig:18s}({src or r.rigidity_confidence:4s}) "
            f"domΔ={_fmt(getattr(r, 'dom_dist_median', None))} domΔcv={_fmt(getattr(r, 'dom_dist_cv', None), 2)} "
            f"n={r.n_structures:>3d} {'human' if r.is_human else 'other':5s} [{tags}] {r.sequence}")
    return "\n".join(lines)


def _fmt(v, nd: int = 1) -> str:
    return "  -  " if v is None or pd.isna(v) else f"{v:.{nd}f}"
