"""Rigid-linker design system: conditional retrieval + evidence calibration + external
structure verification, with linkerdb as the spec, the source and the benchmark.

Why this shape (see scripts/discriminator_experiment.py and design.md section 17): we
tested whether a linker's sequence predicts the database's superposition rigidity under a
homology-grouped split. It does not (AUC ~0.52, chance), and the "rigid" label itself falls
with the number of structures behind it (0.77 at n=3 -> 0.47 at n>20) -- it is partly a
conformational-sampling artifact. So this module does NOT train a sequence->rigidity
classifier. Instead the trustworthy rigidity signal is:

* the fold you design in -- a continuous alpha-helix or beta-strand across the junction is
  a rigid element by construction (``designed_fold``); and
* cross-conformer EVIDENCE from the DB -- a coil called rigid after being seen in many
  structures (ideally different functional states) is trustworthy; a low-n "rigid" is not.

``propose`` turns a design spec into a ranked, evidence-tagged shortlist of natural rigid
linkers (optionally diversified into novel but conservative variants). The real per-construct
call is deferred to external structure prediction; ``verify`` scores what comes back with the
same geometry the DB uses. This narrows and ranks candidates; it does not predict that a
linker will work in a particular fusion -- model and test the shortlist.
"""

from __future__ import annotations

import pandas as pd

from . import developability, superpose
from .artificial import EXTENDED_COIL, HELIX_RISE
from .retrieve import query_canonical

RIGID_LABELS = {"rigid", "rigid_helical", "strand_spanning"}
SPANNING = {"helix_spanning", "strand_spanning", "sse_abutting"}
MIN_STRONG_N = 10   # a coil "rigid" below this #structures may be under-sampled (experiment)

# trust tiers, best first; a candidate's rigidity is only as good as its evidence
TRUST_RANK = {"designed_fold": 0, "evidence_strong": 1, "evidence_weak": 2, "proxy_only": 3}

# conservative, within-class substitutions used for diversification (both directions)
_CONSERVATIVE = {"S": "T", "T": "S", "D": "E", "E": "D", "K": "R", "R": "K",
                 "I": "L", "L": "I", "V": "I", "F": "Y", "Y": "F", "Q": "N", "N": "Q"}


def rigidity_trust(row) -> str | None:
    """Evidence tier for a rigid linker, or None if the row is not rigid.

    ``designed_fold``   continuous SSE across the junction -- rigid by construction.
    ``evidence_strong`` coil, superposition-backed, seen in >= MIN_STRONG_N structures.
    ``evidence_weak``   coil, superposition-backed but few structures (possible under-sampling).
    ``proxy_only``      rigid only by the cheap proxy -- no cross-conformer test at all.
    """
    rig = row.get("rigidity_final") or row.get("rigidity")
    if rig not in RIGID_LABELS:
        return None
    if row.get("linker_class") in SPANNING:
        return "designed_fold"
    if row.get("rigidity_source") == "superposition":
        return "evidence_strong" if (row.get("n_structures") or 0) >= MIN_STRONG_N else "evidence_weak"
    return "proxy_only"


def evidence_tags(row, trust: str) -> list[str]:
    tags = []
    cls = row.get("linker_class")
    if trust == "designed_fold":
        tags.append("continuous-helix rigid rod" if cls == "helix_spanning"
                    else "continuous-strand" if cls == "strand_spanning" else "abutting-SSE")
    elif trust == "evidence_strong":
        tags.append(f"rigid across {int(row.get('n_structures', 0))} structures")
    elif trust == "evidence_weak":
        tags.append(f"UNDER-SAMPLED: rigid but only {int(row.get('n_structures', 0))} structures "
                    f"(low-n rigid is unreliable)")
    else:
        tags.append("proxy-only: no cross-conformer evidence")
    raw = row.get("liability_tags")
    if raw is not None and len(raw):
        tags.extend(list(raw))
    return tags


def reach(length: int, kind: str) -> float:
    """Modelled span (Å): a helix is a rigid rod (HELIX_RISE/res); a coil's upper bound is
    the fully-extended length (EXTENDED_COIL/res)."""
    return (HELIX_RISE if kind == "helix" else EXTENDED_COIL) * length


def _no_new_liability(orig: str, variant: str) -> bool:
    return developability.liability_score(developability.scan(variant)) <= \
        developability.liability_score(developability.scan(orig))


def diversify(seq: str, k: int = 3) -> list[str]:
    """Up to ``k`` conservative single-substitution variants of ``seq`` that introduce no new
    developability liability. Novel sequences that stay inside natural rigid-linker space."""
    seq = (seq or "").upper()
    out, seen = [], {seq}
    for i, aa in enumerate(seq):
        if aa not in _CONSERVATIVE:
            continue
        cand = seq[:i] + _CONSERVATIVE[aa] + seq[i + 1:]
        if cand in seen or not _no_new_liability(seq, cand):
            continue
        seen.add(cand)
        out.append(cand)
        if len(out) >= k:
            break
    return out


def propose(canonical: pd.DataFrame, distance: float | None = None, distance_tol: float = 5.0,
            max_length: int | None = None, min_length: int | None = None, human_only: bool = False,
            no_glyc: bool = False, no_cys: bool = False, max_liability: int | None = None,
            min_trust: str = "evidence_weak", diversify_top: int = 0, top: int = 20) -> pd.DataFrame:
    """A ranked, evidence-tagged shortlist of natural rigid linkers for a design spec.

    ``min_trust`` drops candidates below that evidence tier (default: keep everything at or
    above ``evidence_weak``, i.e. exclude proxy-only). ``diversify_top`` emits conservative
    novel variants for the top N hits.
    """
    hits = query_canonical(canonical, distance=distance, distance_tol=distance_tol, rigidity="rigid",
                           human_only=human_only, min_length=min_length, max_length=max_length,
                           no_glyc=no_glyc, no_cys=no_cys, max_liability=max_liability,
                           top=max(top * 3, 60))
    cutoff = TRUST_RANK.get(min_trust, 2)
    rows = []
    for _, r in hits.iterrows():
        trust = rigidity_trust(r)
        if trust is None or TRUST_RANK[trust] > cutoff:
            continue
        cls = r.get("linker_class")
        coil_seq = r.get("coil_sequence")
        # a coil with its own segment is a splice-in linker; a spanning class is a "join the
        # domains with a continuous SSE" fold motif (the stored sequence is flanking context)
        if cls in SPANNING or not (isinstance(coil_seq, str) and coil_seq):
            use, show_seq = "fold-motif", r.get("sequence")
            kind, length, reach_A = ("helix" if cls == "helix_spanning" else "strand"), 0, None
        else:
            use, show_seq, kind = "splice", coil_seq, "coil"
            length = int(r.get("length_median") or len(coil_seq))
            reach_A = round(reach(length, kind), 1)
        rows.append({
            "candidate_id": r.get("canonical_id"), "sequence": show_seq, "use": use, "origin": "natural",
            "source": r.get("canonical_id"), "linker_class": cls, "length": length,
            "rigidity": r.get("rigidity_final") or r.get("rigidity"), "trust": trust,
            "reach_A": reach_A, "dom_dist": r.get("dom_dist_median"),
            "n_structures": int(r.get("n_structures") or 0), "is_human": bool(r.get("is_human")),
            "liability_score": int(r.get("liability_score") or 0),
            "evidence": evidence_tags(r, trust),
        })
    ranked = sorted(rows, key=lambda d: (TRUST_RANK[d["trust"]], -d["n_structures"], d["liability_score"]))

    out = ranked[:top]
    if diversify_top:
        extra = []
        # only splice-in coil linkers are sensibly mutated; a fold motif's string is context
        for base in [b for b in ranked[:diversify_top] if b["use"] == "splice"]:
            for v in diversify(base["sequence"]):
                nv = dict(base)
                nv.update(candidate_id=f"{base['candidate_id']}~v", sequence=v, origin="variant",
                          evidence=base["evidence"] + [f"conservative variant of {base['source']}"])
                extra.append(nv)
        out = out + extra
    return pd.DataFrame(out)


def format_proposals(df: pd.DataFrame) -> str:
    if df.empty:
        return "no rigid linkers match the spec"
    lines = []
    for r in df.itertuples():
        ev = "; ".join(r.evidence) if isinstance(r.evidence, list) else str(r.evidence)
        span = f"reach~{r.reach_A:>5.1f}A" if r.reach_A is not None \
            else (f"domΔ~{r.dom_dist:>5.1f}A" if pd.notna(r.dom_dist) else "domΔ  -  ")
        lines.append(f"{r.candidate_id:24s} {r.origin:7s} {r.use:10s} len={r.length:>3d} {r.linker_class:14s} "
                     f"{r.rigidity:14s} trust={r.trust:14s} {span} n={r.n_structures:>3d} "
                     f"{'human' if r.is_human else 'other':5s} {r.sequence}")
        lines.append(f"    evidence: {ev}")
    lines.append("\nuse=splice: a coil linker to graft in; use=fold-motif: join the domains with a continuous")
    lines.append("SSE (the string is flanking context, not a splice-in sequence). Ranking is by evidence tier,")
    lines.append("not a predicted success probability -- the rigidity of a novel linker in YOUR fusion is set by")
    lines.append("the two domains and the fold; verify externally (design.FOLD_PROTOCOL).")
    return "\n".join(lines)


# ---- external structure verification -------------------------------------------------
# linkerdb cannot fold proteins. The per-construct rigidity call must come from an external
# structure predictor; this module scores what it returns, with the DB's own geometry.

FOLD_PROTOCOL = """\
External verification protocol (run outside linkerdb; GPU needed):
 1. Build the fusion sequence: domain_N + candidate_linker + domain_C.
 2. Predict an ENSEMBLE of structures, not one model -- rigidity is a cross-conformer
    property (see design.md section 17). Options:
      * AlphaFold-multimer / AlphaFold3 or ESMFold with several seeds / MSA subsamplings;
      * RFdiffusion inpainting of the linker between the two fixed domains, many samples.
 3. For each model, extract Cα of the two domains and feed them to design.ensemble_rigidity,
    which applies the SAME superposition measure (c_rmsd / hinge) the database is built on.
 4. Feasibility of a single model: design.spans() checks the linker physically bridges the
    intended domain separation. Rigidity: ensemble_rigidity over the models.
"""


def spans(linker_len: int, kind: str, end_to_end_A: float) -> bool:
    """Is an observed/target end-to-end distance within what a linker of this length can span?
    (helix: rigid-rod length; coil: fully-extended upper bound.)"""
    return end_to_end_A <= reach(linker_len, kind) + 1e-9


def ensemble_rigidity(models: list[tuple[str, dict, dict]]) -> dict:
    """Rigidity of a predicted ensemble, using the DB's superposition measure.

    ``models``: list of (label, n_domain_ca, c_domain_ca); each ``*_ca`` is
    ``{residue_index: xyz}`` with indices shared across models (reference first). Returns the
    superposition ``max_c_rmsd`` / ``max_hinge`` and the ``rigid`` / ``flexible`` /
    ``intermediate`` / ``undetermined`` call -- identical thresholds to superpose.classify.
    """
    mob = superpose.mobility(models)
    mob["rigidity"] = superpose.classify(mob["max_c_rmsd"], mob["max_hinge"])
    return mob
