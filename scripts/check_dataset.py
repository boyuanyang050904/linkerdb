"""Sanity checks on a built dataset directory (default: results/full).

Usage: python scripts/check_dataset.py [OUT_DIR]

Every check prints PASS/FAIL; the script exits non-zero if any invariant is violated.
"""

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

out = Path(sys.argv[1] if len(sys.argv) > 1 else "results/full")
final = pd.read_parquet(out / "linkers.parquet")
allk = pd.read_parquet(out / "linkers_all.parquet")
chains = pd.read_parquet(out / "chains.parquet")
info = json.loads((out / "run_info.json").read_text())
max_res = info["rules"]["max_resolution"]
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


excluding = {"long_gap", "unp_unmapped", "unp_discontinuous", "engineered_linker", "geometry_violation",
             "overlapping_segments", "no_uniprot", "chimera", "not_protein"}
check("all final linkers pass QC", final["qc_pass"].all())
check("no excluding flag in final set", not final["qc_flags"].map(lambda f: bool(set(f) & excluding)).any())
check(f"resolution <= {max_res} Å", (final["resolution"] <= max_res).all(), f"max {final['resolution'].max()}")
check("no NMR entries", not final["method"].str.contains("NMR").any())
check("ss_length == len(ss_sequence)", (final["ss_length"] == final["ss_sequence"].str.len()).all())
check("ss_length == len(ss_dssp)", (final["ss_length"] == final["ss_dssp"].str.len()).all())
check("gap inside structural linker", ((final["ss_start"] <= final["gap_start"]) & (final["gap_end"] <= final["ss_end"])
                                       | (final["gap_length"] == 0)).all())
lk = final[final["ss_length"] > 0]
check("structural linkers contain no regular SSE outside the CATH gap",
      all(not any(c in "HEI" for i, c in enumerate(r.ss_dssp)
                  if not (r.gap_start <= r.ss_start + i <= r.gap_end))
          for r in lk.itertuples()))
e2e = lk.dropna(subset=["end_to_end"])
check("end-to-end distance physically possible",
      (e2e["end_to_end"] <= 3.8 * (e2e["ss_length"] - 1) + 0.5).all())
span = final.dropna(subset=["ca_span"])
check("anchor span physically possible", (span["ca_span"] <= 3.8 * (span["ss_length"] + 1) + 0.5).all())
check("one representative chain per 40% cluster",
      final.groupby("cluster")[["pdb_id", "auth_asym_id"]].nunique().max().max() == 1)
check("linker ids unique", final["linker_id"].is_unique)
# no residual same-UniProt overlaps across chains
dup = 0
for acc, g in final.dropna(subset=["uniprot"]).groupby("uniprot"):
    rows = list(g.itertuples())
    for i, a in enumerate(rows):
        for b in rows[i + 1:]:
            if (a.pdb_id, a.auth_asym_id) != (b.pdb_id, b.auth_asym_id) and a.unp_lo <= b.unp_hi and b.unp_lo <= a.unp_hi:
                dup += 1
check("no cross-chain overlap of UniProt linker regions", dup == 0, f"{dup} pairs")
check("excluded chimeras never reach the final set",
      not set(final["pdb_id"] + final["auth_asym_id"]) & set(
          chains.loc[chains["status"].str.startswith("excluded"), "pdb_id"]
          + chains.loc[chains["status"].str.startswith("excluded"), "auth_asym_id"].fillna("")))

print()
print("chains:", dict(chains["status"].value_counts()))
print("final linkers:", len(final), "from", final.groupby(["pdb_id", "auth_asym_id"]).ngroups, "chains")
print("dedup dropped:", Counter(str(x).split(":")[0] for x in allk["dedup_dropped"].dropna()))
print("flags in final:", Counter(f for fl in final["qc_flags"] for f in fl).most_common())
print("junction_ss:", dict(final["junction_ss"].value_counts()))
print("junction_type:", dict(final["junction_type"].value_counts()))
nz = final[final["ss_length"] > 0]["ss_length"]
print("structural linker length (non-empty): mean %.1f median %d q90 %d" % (nz.mean(), nz.median(), nz.quantile(.9)))
print("gap==0: %.1f%%  ss==0: %.1f%%" % (100 * np.mean(final["gap_length"] == 0), 100 * np.mean(final["ss_length"] == 0)))
sys.exit(1 if fails else 0)
