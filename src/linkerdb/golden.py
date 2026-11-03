"""Golden set: a small, hand-curated reference of well-known multi-domain proteins
whose inter-domain rigidity is textbook, used to sanity-check (and eventually
calibrate) the database's rigidity call.

Honesty note: this is a *seed* set (~a dozen entries) assembled from established
structural-biology knowledge, not a literature-curated corpus with citations, and
it carries no hand-verified residue boundaries. Every identifier is verified to be
present in the database, but ``expected_rigidity`` is a qualitative judgement. It is
enough to expose systematic errors (it already shows that domain-centroid distance
misses hinge *rotation*); it is NOT enough to tune thresholds rigorously. Expand it
with expert-curated entries before using it to set parameters.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

GOLDEN_CSV = Path(__file__).resolve().parents[2] / "golden" / "golden_set.csv"

# our fine-grained rigidity labels -> the coarse rigid/flexible axis the golden set uses
COARSE = {
    "rigid": "rigid", "rigid_helical": "rigid", "strand_spanning": "rigid",
    "flexible": "flexible", "flexible_disordered": "flexible",
    "intermediate": "intermediate", "undetermined": "undetermined",
}


def coarse(label: str) -> str:
    return COARSE.get(label, "undetermined")


def load_golden(path: Path | str = GOLDEN_CSV) -> pd.DataFrame:
    return pd.read_csv(path)


def score_golden(golden: pd.DataFrame, canonical: pd.DataFrame) -> pd.DataFrame:
    """For each golden protein, match its best-supported canonical linker and compare
    our rigidity (coarsened) to the expected one. One row per golden entry."""
    rows = []
    for e in golden.itertuples(index=False):
        m = canonical[canonical["uniprot"] == e.uniprot]
        if "region_lo" in golden.columns and pd.notna(getattr(e, "region_lo", None)):
            m = m[(m["unp_lo"] <= e.region_hi) & (m["unp_hi"] >= e.region_lo)]
        if m.empty:
            rows.append({"protein": e.protein, "uniprot": e.uniprot, "expected": e.expected_rigidity,
                         "confidence": e.confidence, "found": False, "ours": None, "rigidity": None,
                         "match": None})
            continue
        top = m.sort_values("n_structures", ascending=False).iloc[0]
        label = top["rigidity_final"] if "rigidity_final" in m.columns and pd.notna(top.get("rigidity_final")) \
            else top["rigidity"]
        ours = coarse(label)
        match = None if ours == "intermediate" else (ours == e.expected_rigidity)
        rows.append({
            "protein": e.protein, "uniprot": e.uniprot, "expected": e.expected_rigidity,
            "confidence": e.confidence, "found": True, "ours": ours, "rigidity": label,
            "canonical_id": top["canonical_id"], "n_structures": int(top["n_structures"]),
            "rigidity_confidence": top["rigidity_confidence"],
            "dom_dist_cv": top.get("dom_dist_cv"), "elbow_angle_std": top.get("elbow_angle_std"),
            "elbow_angle_range": top.get("elbow_angle_range"),
            "match": match,
        })
    return pd.DataFrame(rows)


def summary(scored: pd.DataFrame) -> str:
    found = scored[scored["found"]]
    decided = found[found["match"].notna()]
    n_match = int(decided["match"].sum())
    lines = [f"golden set: {len(scored)} entries, {len(found)} found in DB, "
             f"{len(found) - len(decided)} scored 'intermediate' (no credit either way)",
             f"decided: {n_match}/{len(decided)} agree with expected rigidity"
             + (f" ({n_match / len(decided):.0%})" if len(decided) else "")]
    for r in scored.itertuples(index=False):
        if not r.found:
            lines.append(f"  {r.protein:28s} {r.uniprot}  NOT FOUND in DB")
            continue
        mark = {True: "OK ", False: "XX ", None: "~~ "}[r.match]
        lines.append(f"  {mark}{r.protein:28s} expected={r.expected:8s} ours={r.rigidity:16s} "
                     f"(n={r.n_structures}, {r.rigidity_confidence})")
    return "\n".join(lines)
