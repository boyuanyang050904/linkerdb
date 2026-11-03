"""Summary statistics and a human-readable view of single linkers for manual checks."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

AA20 = "ACDEFGHIKLMNPQRSTVWY"


def summarize_canonical(canonical: pd.DataFrame, instances: pd.DataFrame | None = None,
                        title: str = "canonical linkers") -> str:
    lines = [f"===== {title} ====="]
    if instances is not None and not instances.empty:
        lines.append(f"instances: {len(instances)} junctions, {int(instances['qc_pass'].sum())} qc_pass "
                     f"from {instances.groupby(['pdb_id', 'auth_asym_id']).ngroups} chains")
    if canonical.empty:
        lines.append("no canonical linkers")
        return "\n".join(lines)
    c = canonical
    grouped = c[c["grouped"]]
    lines.append(f"canonical linkers: {len(c)}  (grouped by UniProt: {len(grouped)}, ungrouped: {len(c) - len(grouped)})")
    lines.append(f"distinct UniProt: {c['uniprot'].nunique()}   human: {int(c['is_human'].sum())} "
                 f"({c['is_human'].mean():.0%})")
    multi = c[c["n_structures"] >= 3]
    lines.append(f"with >=3 structures (rigidity assessable): {len(multi)} ({len(multi) / len(c):.0%}); "
                 f"max structures for one linker: {int(c['n_structures'].max())}")
    lines.append("rigidity: " + ", ".join(f"{k}={v}" for k, v in c["rigidity"].value_counts().items()))
    lines.append("rigidity confidence: " + ", ".join(f"{k}={v}" for k, v in
                                                      c["rigidity_confidence"].value_counts().items()))
    lines.append("length (median per linker): " + _q(c["length_median"]))
    hi = c[c["end_to_end_cv"].notna()]
    if len(hi):
        lines.append(f"end-to-end CV (flexibility proxy, n>=2 structures): median {hi['end_to_end_cv'].median():.2f} "
                     f"over {len(hi)} linkers")
    liab = {k: int((c[k] > 0).sum()) for k in ("n_glyc", "deamid_ng", "isomer_dg", "n_cys", "furin") if k in c}
    lines.append("linkers carrying a liability: " + ", ".join(f"{k}={v} ({v / len(c):.0%})" for k, v in liab.items()))
    clean = c[(c.get("n_glyc", 0) == 0) & (c.get("deamid_ng", 0) == 0) & (c.get("isomer_dg", 0) == 0) &
              (c.get("furin", 0) == 0)]
    lines.append(f"liability-free linkers: {len(clean)} ({len(clean) / len(c):.0%})")
    return "\n".join(lines)


def _q(s: pd.Series) -> str:
    if s.empty:
        return "n/a"
    q = s.quantile([0, 0.25, 0.5, 0.75, 0.95, 1.0]).tolist()
    return "min {:.0f} | q25 {:.0f} | median {:.0f} | q75 {:.0f} | q95 {:.0f} | max {:.0f} | mean {:.1f}".format(
        *q, s.mean())


def composition(seqs) -> dict[str, float]:
    c = Counter("".join(seqs))
    tot = sum(c[a] for a in AA20)
    return {a: c[a] / tot for a in AA20} if tot else {}


def summarize(chains: pd.DataFrame, linkers: pd.DataFrame, background: pd.Series | None = None,
              title: str = "summary") -> str:
    lines = [f"===== {title} ====="]
    if not chains.empty:
        lines.append(f"entries: {chains['pdb_id'].nunique()}   chains: {len(chains)}")
        lines.append("chain status: " + ", ".join(f"{k}={v}" for k, v in chains["status"].value_counts().items()))
        if "numbering_check" in chains:
            nc = chains["numbering_check"].dropna().value_counts()
            lines.append("numbering check: " + ", ".join(f"{k}={v}" for k, v in nc.items()))
        if "method" in chains:
            m = chains.dropna(subset=["method"])["method"].value_counts()
            lines.append("methods: " + ", ".join(f"{k}={v}" for k, v in m.items()))
        if "dssp_vs_annotation" in chains:
            agr = chains["dssp_vs_annotation"].dropna()
            if len(agr):
                lines.append(f"DSSP vs mmCIF SS annotation (3-state): median {agr.median():.3f}, "
                             f"min {agr.min():.3f} over {len(agr)} chains")
    if linkers.empty:
        lines.append("no linkers")
        return "\n".join(lines)
    ok = linkers[linkers["qc_pass"]]
    lines.append(f"junctions/linkers: {len(linkers)}   qc_pass: {len(ok)}")
    flags = Counter(f for fl in linkers["qc_flags"] for f in fl)
    lines.append("qc flags: " + ", ".join(f"{k}={v}" for k, v in flags.most_common()))
    lines.append("junction type: " + ", ".join(f"{k}={v}" for k, v in ok["junction_type"].value_counts().items()))
    lines.append("junction SS: " + ", ".join(f"{k}={v}" for k, v in ok["junction_ss"].value_counts().items()))
    lines.append(f"gap length (CATH, qc_pass): {_q(ok['gap_length'])};  zero-length {np.mean(ok['gap_length'] == 0):.1%}")
    lines.append(f"structural linker length (qc_pass): {_q(ok['ss_length'])};  zero-length {np.mean(ok['ss_length'] == 0):.1%}")
    nz = ok[ok["ss_length"] > 0]
    if len(nz):
        lines.append(f"  non-empty structural linkers: {len(nz)}; with missing residues {np.mean(nz['n_missing'] > 0):.1%}; "
                     f"fully disordered {np.mean(nz['frac_observed'] == 0):.1%}; capped {np.mean(nz['qc_flags'].map(lambda f: 'capped_n' in f or 'capped_c' in f)):.1%}")
        lines.append(f"  with chain break {np.mean(nz['n_breaks'] > 0):.1%}; with non-standard residue "
                     f"{np.mean(nz['n_nonstandard'] > 0):.1%} ({Counter(x for l in nz['nonstandard'] for x in l).most_common(5)})")
        comp = composition(nz["ss_sequence"])
        if background is not None and len(background):
            bg = composition(background)
            enr = sorted(((comp[a] / bg[a] if bg.get(a) else np.nan, a) for a in AA20), reverse=True)
            lines.append("  linker composition enrichment vs whole chains: " +
                         " ".join(f"{a}:{r:.2f}" for r, a in enr))
    pairs = Counter(zip(ok["cath_n"], ok["cath_c"])).most_common(5)
    lines.append("top CATH pairs: " + "; ".join(f"{a}->{b} ({n})" for (a, b), n in pairs))
    return "\n".join(lines)


def show_linker(chain, rec: dict, ss: str, pad: int = 15) -> str:
    """Text view around one junction: sequence, domain track, DSSP, observed track."""
    lo = max(1, min(rec["ss_start"], rec["gap_start"]) - pad)
    hi = min(chain.n, max(rec["ss_end"], rec["gap_end"]) + pad)
    seq = "".join(chain.res(i).one for i in range(lo, hi + 1))
    dom = []
    for i in range(lo, hi + 1):
        if i < rec["gap_start"]:
            dom.append("N")
        elif i > rec["gap_end"]:
            dom.append("C")
        else:
            dom.append("-")
    lk = "".join("L" if rec["ss_start"] <= i <= rec["ss_end"] else "." for i in range(lo, hi + 1))
    obs = "".join("|" if chain.res(i).observed else "x" for i in range(lo, hi + 1))
    return "\n".join([
        f"{rec['linker_id']}  {rec['domain_n']}({rec['cath_n']}) -> {rec['domain_c']}({rec['cath_c']})  "
        f"type={rec['junction_type']} gap={rec['gap_length']} ss_len={rec['ss_length']} "
        f"junction_ss={rec['junction_ss']} flags={','.join(rec['qc_flags'])}",
        f"  SEQRES {lo}-{hi} (auth {chain.res(lo).auth}-{chain.res(hi).auth})",
        "  seq  " + seq,
        "  CATH " + "".join(dom),
        "  DSSP " + ss[lo - 1:hi],
        "  link " + lk,
        "  obs  " + obs,
        f"  linker: {rec['ss_sequence'] or '(empty)'}  SEQRES {rec['ss_start']}-{rec['ss_end']}"
        f"  auth {rec['ss_start_auth'] or '-'}-{rec['ss_end_auth'] or '-'}",
    ])
