"""Sequence-based developability liabilities of a linker.

These are *sequence motifs* that are known chemical/biological liability sites,
not predictions that a given linker will fail: whether a motif is actually a
problem depends on solvent exposure, local structure and the expression system.
A linker, being exposed and flexible, is exactly where these sites are most
likely to be reactive, which is why they are worth flagging for a designer.

Flags (counts unless noted):

* ``n_glyc``       N-linked glycosylation sequon  N-X-[S/T], X != P
* ``deamid_ng``    strong deamidation hotspot     N-G
* ``deamid_soft``  softer deamidation             N-[S/T/H/A]
* ``isomer_dg``    aspartate isomerization        D-G
* ``frag_ds``      Asp fragmentation / cleavage   D-[S/P]
* ``n_met``, ``n_trp``   oxidation-prone residues
* ``n_cys``        cysteines (free-thiol / scrambling risk; SS-bond state unknown from sequence)
* ``furin``        furin-like cleavage            R-X-[R/K]-R
* ``n_lys``        extra lysines (chemical-conjugation / charge heterogeneity handle)
"""

from __future__ import annotations

import re

_PATTERNS = {
    "n_glyc": re.compile(r"N[^P][ST]"),
    "deamid_ng": re.compile(r"NG"),
    "deamid_soft": re.compile(r"N[STHA]"),
    "isomer_dg": re.compile(r"DG"),
    "frag_ds": re.compile(r"D[SP]"),
    "furin": re.compile(r"R.[RK]R"),
}
_RESIDUE_COUNTS = {"n_met": "M", "n_trp": "W", "n_cys": "C", "n_lys": "K"}

LIABILITY_KEYS = tuple(_PATTERNS) + tuple(_RESIDUE_COUNTS)


def scan(seq: str) -> dict[str, int]:
    """Count liability motifs in ``seq`` (one-letter, upper-case; 'X'/unknowns ignored).

    Overlapping matches are counted (e.g. ``NN`` in ``NNS`` -> one sequon at the
    first N), using look-ahead so adjacent sites are not missed.
    """
    s = (seq or "").upper()
    out = {k: len(re.findall(f"(?={p.pattern})", s)) for k, p in _PATTERNS.items()}
    out.update({k: s.count(aa) for k, aa in _RESIDUE_COUNTS.items()})
    return out


def summary_flags(counts: dict[str, int]) -> list[str]:
    """Human-readable liability tags for the motifs that are present."""
    tags = []
    if counts.get("n_glyc"):
        tags.append("N-glycosylation")
    if counts.get("deamid_ng"):
        tags.append("deamidation(NG)")
    if counts.get("isomer_dg"):
        tags.append("isomerization(DG)")
    if counts.get("furin"):
        tags.append("furin-site")
    if counts.get("n_cys"):
        tags.append("cysteine")
    if counts.get("n_met", 0) + counts.get("n_trp", 0) >= 2:
        tags.append("oxidation-prone")
    return tags


def liability_score(counts: dict[str, int]) -> int:
    """A crude additive count of hard liabilities, for ranking only (not a risk model)."""
    return (3 * counts.get("n_glyc", 0) + 2 * counts.get("deamid_ng", 0) + 2 * counts.get("isomer_dg", 0)
            + 3 * counts.get("furin", 0) + counts.get("n_cys", 0) + counts.get("frag_ds", 0))
