"""Retrospective validation against published linker-swap experiments.

The golden set (``golden.py``) asks a different question from this module. There we
check whether our *rigidity label* for a natural multi-domain protein matches its
textbook conformational behaviour. Here we check whether the library's **mechanistic
models** -- the ones a fusion-protein designer actually leans on -- line up with what
happened when someone swapped the linker in a real engineered construct and measured
activity, expression or oligomeric state.

The library makes exactly two kinds of falsifiable claim that a linker swap can test:

* **reach / geometry** -- how far a linker holds two domains apart. A helix extends
  ``HELIX_RISE`` (~1.5 A) per residue along its axis; a flexible coil has no fixed
  length, with a fully-extended *upper bound* of ``EXTENDED_COIL`` (~3.5 A) per residue
  (both from ``artificial.py``). From these the library predicts, e.g., the minimum
  flexible-linker length needed to span a given gap.
* **rigidity** -- whether a linker holds a fixed relative geometry (helix / EAAAK-like)
  or samples a broad distribution (coil / GS-like).

A linker swap often turns on things the library says nothing about -- expression yield,
proteolysis, the specific biology of the fused partners. Those rows are scored
``no_geometry_prediction`` and excluded from the agreement rate, the same way the
golden set gives no credit for ``intermediate``. Pretending the library can adjudicate
them would be dishonest; showing where it is silent is part of the validation.

Non-circularity: ``golden/linker_swap_set.csv`` stores *only the experimental outcome*.
The library's prediction for each row is computed here from the library's own constants
and class convention, then compared. Changing ``HELIX_RISE`` / ``EXTENDED_COIL`` moves
the score, which is the coupling we want.

Honesty note (mirrors ``golden.py``): the bundled corpus is a *seed* assembled from
established protein-engineering knowledge. The ``reference`` strings are recollections
flagged for verification, not checked citations, and there are no hand-verified residue
boundaries or digitised measurements. It is enough to show whether the models point the
right way on textbook cases; it is NOT a citation-curated benchmark. Expand it with
primary literature (and the measured distances / lengths) before drawing quantitative
conclusions. See design.md section 16.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from .artificial import EXTENDED_COIL, HELIX_RISE

SWAP_CSV = Path(__file__).resolve().parents[2] / "golden" / "linker_swap_set.csv"

# what rigidity class the library assigns a linker *type*, and a rigid>flexible rank
TYPE_CLASS = {"helix": "rigid_helical", "strand": "strand_spanning", "coil": "flexible", "other": "unknown"}
TYPE_RANK = {"helix": 2, "strand": 2, "coil": 0, "other": None}  # None = library has no call


def load_swaps(path: Path | str = SWAP_CSV) -> pd.DataFrame:
    return pd.read_csv(path)


def _num(v) -> float | None:
    """Parse a possibly-blank cell to float, else None."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def predict(row) -> tuple[str, bool | None]:
    """Compute the library's prediction for one swap row and whether it agrees with the
    experiment. Returns ``(prediction_text, agree)`` where ``agree`` is None when the
    library makes no prediction on this row (out of scope, or an undistinguishable tie).
    """
    kind = row["expected_kind"]
    la_type = row.get("la_type")
    la_len = _num(row.get("la_len"))
    target = _num(row.get("span_target_A"))
    tol = _num(row.get("tol"))

    if kind == "min_len_for_span":
        # fewest coil residues that could bridge the gap (fully-extended lower bound)
        pred_len = math.ceil(target / EXTENDED_COIL)
        obs = _num(row["expected"])
        agree = abs(pred_len - obs) <= (tol if tol is not None else 0)
        return (f"coil needs >={pred_len} res to span {target:.0f} A (3.5 A/res); "
                f"observed threshold ~{obs:.0f}"), agree

    if kind in ("span_insufficient", "span_sufficient"):
        reach_ub = EXTENDED_COIL * la_len  # coil upper bound; helix would use HELIX_RISE
        if la_type == "helix":
            reach_ub = HELIX_RISE * la_len
        can_span = reach_ub >= target
        want_span = kind == "span_sufficient"
        agree = can_span == want_span
        return (f"{la_type} len {la_len:.0f} max reach {reach_ub:.0f} A "
                f"{'>=' if can_span else '<'} {target:.0f} A -> "
                f"{'can' if can_span else 'cannot'} span"), agree

    if kind == "length_distance_coupling":
        # rigid rod: length sets distance; flexible coil: it does not
        pred = "controls" if TYPE_CLASS.get(la_type, "").startswith("rigid") or la_type == "strand" else "weak"
        agree = pred == str(row["expected"]).strip()
        return f"{la_type}: length->distance coupling predicted '{pred}'", agree

    if kind == "more_rigid":
        ra, rb = TYPE_RANK.get(la_type), TYPE_RANK.get(row.get("lb_type"))
        if ra is None or rb is None:
            return "library has no rigidity call for one side", None
        pred = "a" if ra > rb else ("b" if rb > ra else "tie")
        if pred == "tie":
            return "both sides predicted equally rigid (tie)", None
        agree = pred == str(row["expected"]).strip()
        return (f"more rigid = {pred} ({row.get('la_name') if pred == 'a' else row.get('lb_name')}, "
                f"{TYPE_CLASS.get(la_type if pred == 'a' else row.get('lb_type'))})"), agree

    if kind == "no_geometry_prediction":
        return "out of scope: library predicts geometry/rigidity, not this observable", None

    return f"unknown expected_kind '{kind}'", None


def score_swaps(swaps: pd.DataFrame) -> pd.DataFrame:
    """One row per corpus entry with the library prediction and the agree verdict."""
    rows = []
    for _, e in swaps.iterrows():
        pred_text, agree = predict(e)
        rows.append({
            "system": e["system"], "axis": e["axis"], "expected_kind": e["expected_kind"],
            "observable": e["observable"], "expected": e["expected"],
            "prediction": pred_text, "agree": agree,
            "scored": agree is not None, "confidence": e.get("confidence"),
            "reference": e.get("reference"),
        })
    return pd.DataFrame(rows)


def summary(scored: pd.DataFrame) -> str:
    decided = scored[scored["scored"]]
    coverage = scored[~scored["scored"]]
    n_match = int(decided["agree"].sum()) if not decided.empty else 0
    lines = [
        "linker-swap retrospective (SEED corpus -- recollection-based, not citation-curated; "
        "see design.md section 16)",
        f"  {len(scored)} comparisons: {len(decided)} scored on the reach/rigidity axes, "
        f"{len(coverage)} out of scope (geometry cannot adjudicate -> no credit)",
        f"  scored agreement: {n_match}/{len(decided)}"
        + (f" ({n_match / len(decided):.0%})" if len(decided) else ""),
    ]
    for r in scored.itertuples(index=False):
        mark = {True: "OK ", False: "XX ", None: "~~ "}[r.agree]
        tag = f"[{r.axis}]"
        lines.append(f"  {mark}{tag:11s} {r.system}")
        lines.append(f"       exp: {r.observable}")
        lines.append(f"       lib: {r.prediction}")
    lines.append("  OK/XX = scored and agrees/disagrees; ~~ = out of scope or undecidable (not scored)")
    return "\n".join(lines)
