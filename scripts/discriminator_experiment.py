"""Decision-gate experiment for the rigid-linker design system (design.md section 17).

Question: can a linker's SEQUENCE predict the database's superposition rigidity label,
well enough to serve as a learned discriminator for de novo linkers?

Method: logistic regression (numpy only) on sequence-derived features, evaluated with
cross-validation GROUPED BY UniProt accession so that homologous linkers never straddle
the train/test split (a random split leaks homology and inflates the score). Compared
against a length-only model and a leakage upper bound that adds structure-derived columns,
and against a random (leaky) split for contrast.

Result (DB built 2026-11-03): sequence-only AUC ~= 0.52 under the grouped split -- no
better than the 0.71 base rate / chance. Even the structure-derived upper bound reaches
only ~0.61. Crucially, rigid-fraction FALLS with n_structures (0.77 at n=3 -> 0.47 at
n>20): the "rigid" label is substantially a conformational-SAMPLING artifact, not a
sequence property. So the design system does NOT train a sequence->rigidity classifier;
it uses the fold you design in (a continuous helix is a rigid rod) plus cross-conformer
EVIDENCE strength from the DB, and defers the real call to external structure prediction.
This script is the evidence for that choice; re-run it to reproduce the numbers.

    python scripts/discriminator_experiment.py results/db/linkers_canonical.parquet
"""
from __future__ import annotations

import collections
import sys

import numpy as np
import pandas as pd

from linkerdb import developability

AA = "ACDEFGHIKLMNPQRSTVWY"
# Pace & Scholtz (1998) helix propensity, kcal/mol (lower = stronger former); negated so higher = more helical.
PACE = {"A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40, "I": 0.41,
        "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61, "N": 0.65, "T": 0.66,
        "C": 0.68, "D": 0.69, "G": 1.00, "P": 3.16}


def seq_feats(seq: str) -> dict:
    s = "".join(c for c in (seq or "").upper() if c in AA)
    n = max(len(s), 1)
    f = {"loglen": np.log1p(len(s))}
    for name, res in [("f_gly", "G"), ("f_pro", "P"), ("f_ala", "A"), ("f_ser", "S")]:
        f[name] = s.count(res) / n
    f["f_charged"] = sum(s.count(c) for c in "DEKR") / n
    f["f_arom"] = sum(s.count(c) for c in "FWY") / n
    f["f_hphob"] = sum(s.count(c) for c in "AILMFWV") / n
    f["f_gp"] = (s.count("G") + s.count("P")) / n
    f["helix_prop"] = -np.mean([PACE[c] for c in s]) if s else 0.0
    f["net_charge"] = (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / n
    f["liab"] = developability.liability_score(developability.scan(s))
    return f


SEQ_COLS = ["loglen", "f_gly", "f_pro", "f_ala", "f_ser", "f_charged", "f_arom", "f_hphob",
            "f_gp", "helix_prop", "net_charge", "liab"]


def _standardize(Xtr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1
    return (Xtr - mu) / sd, (Xte - mu) / sd


def fit_logreg(X, y, l2=1.0, iters=50):
    """L2-regularized logistic regression by IRLS/Newton (bias unpenalized)."""
    X = np.hstack([np.ones((len(X), 1)), X])
    w = np.zeros(X.shape[1])
    R = np.eye(X.shape[1]) * l2
    R[0, 0] = 0
    for _ in range(iters):
        p = np.clip(1 / (1 + np.exp(-X @ w)), 1e-9, 1 - 1e-9)
        H = X.T * (p * (1 - p)) @ X + R
        g = X.T @ (y - p) - R @ w
        try:
            w += np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
    return w


def _predict(w, X):
    return 1 / (1 + np.exp(-np.hstack([np.ones((len(X), 1)), X]) @ w))


def auc(y, s):
    pos, neg = (y == 1).sum(), (y == 0).sum()
    if pos == 0 or neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1)
    sr, c = collections.defaultdict(float), collections.defaultdict(int)
    _, inv = np.unique(s, return_inverse=True)
    for r, v in zip(ranks, inv):
        sr[v] += r
        c[v] += 1
    adj = np.array([sr[v] / c[v] for v in inv])
    return (adj[y == 1].sum() - pos * (pos + 1) / 2) / (pos * neg)


def _grouped_folds(groups, k=5, seed=0):
    uniq = pd.unique(groups)
    np.random.RandomState(seed).shuffle(uniq)
    a = {g: i % k for i, g in enumerate(uniq)}
    return np.array([a[g] for g in groups])


def _random_folds(n, k=5, seed=0):
    idx = np.random.RandomState(seed).permutation(n)
    f = np.empty(n, int)
    f[idx] = np.arange(n) % k
    return f


def _cv(df, cols, y, foldid, l2=1.0):
    X = df[cols].to_numpy(float)
    oof = np.zeros(len(df))
    for k in np.unique(foldid):
        tr, te = foldid != k, foldid == k
        Xtr, Xte = _standardize(X[tr], X[te])
        oof[te] = _predict(fit_logreg(Xtr, y[tr], l2=l2), Xte)
    return oof


def _report(name, y, oof):
    print(f"  {name:40s} AUC={auc(y, oof):.3f}  acc@.5={(((oof >= .5) == y).mean()):.3f}  "
          f"Brier={np.mean((oof - y) ** 2):.3f}")
    return auc(y, oof)


def main(path="results/db/linkers_canonical.parquet"):
    df = pd.read_parquet(path)
    sp = df[df["superpose_rigidity"].isin(["rigid", "flexible"])].copy().reset_index(drop=True)
    y = (sp["superpose_rigidity"] == "rigid").astype(int).to_numpy()
    base = y.mean()
    print(f"n={len(sp)}  uniprot groups={sp['uniprot'].nunique()}  "
          f"rigid rate={base:.3f} (majority-class acc={max(base, 1 - base):.3f})\n")

    sp = pd.concat([sp, pd.DataFrame([seq_feats(s) for s in sp["sequence"]])], axis=1)
    sp["is_helix_sp"] = (sp["linker_class"] == "helix_spanning").astype(float)
    sp["is_strand_sp"] = (sp["linker_class"] == "strand_spanning").astype(float)
    sp["dom_dist_median"] = sp["dom_dist_median"].fillna(sp["dom_dist_median"].median())
    struct_cols = SEQ_COLS + ["frac_helix", "frac_disordered", "is_helix_sp", "is_strand_sp", "dom_dist_median"]

    gf, rf = _grouped_folds(sp["uniprot"].to_numpy()), _random_folds(len(sp))
    print("GROUPED-by-UniProt 5-fold CV (honest, no homology leakage):")
    _report("length only", y, _cv(sp, ["loglen"], y, gf))
    _report("sequence-only (comp + helix + len + liab)", y, _cv(sp, SEQ_COLS, y, gf))
    _report("+ structure-derived (LEAKAGE upper bound)", y, _cv(sp, struct_cols, y, gf))
    print("\nRANDOM 5-fold CV (homology leaks across folds, for contrast):")
    _report("sequence-only (random split)", y, _cv(sp, SEQ_COLS, y, rf))

    print("\nThe sampling confound -- rigid fraction vs #structures behind the call:")
    sp["nb"] = pd.cut(sp["n_structures"], [2, 3, 5, 10, 20, 10000],
                      labels=["3", "4-5", "6-10", "11-20", ">20"])
    for b, g in sp.groupby("nb", observed=True):
        print(f"    n_structures {b:>5s}: rigid={g['rig' if 'rig' in g else 'superpose_rigidity'].eq('rigid').mean():.2f}  (N={len(g)})")
    print("\n=> sequence does not carry the rigidity label; low-n 'rigid' is under-sampled.")
    print("   design.py therefore uses designed-fold geometry + DB evidence strength, not a learned classifier.")


if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
