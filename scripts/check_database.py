"""Invariant checks on a built linker database (default: results/db).

Usage: python scripts/check_database.py [DB_DIR]
"""

import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

d = Path(sys.argv[1] if len(sys.argv) > 1 else "results/db")
inst = pd.read_parquet(d / "linkers_instances.parquet")
canon = pd.read_parquet(d / "linkers_canonical.parquet")
chains = pd.read_parquet(d / "chains.parquet")
info = json.loads((d / "run_info.json").read_text())
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


ok = inst[inst["qc_pass"]]
check("instances: qc_pass resolution <= 2.5", (ok["resolution"] <= 2.5).all(), f"max {ok['resolution'].max()}")
check("instances: no NMR", not ok["method"].str.contains("NMR").any())
check("canonical ids unique", canon["canonical_id"].is_unique)
check("linker_class vocabulary",
      set(canon["linker_class"]) <= {"coil", "helix_spanning", "strand_spanning", "sse_abutting"})
check("coil linkers carry a coil sequence",
      (canon.loc[canon["linker_class"] == "coil", "coil_sequence"].str.len() > 0).all())
check("grouped linkers have a UniProt interval",
      canon[canon["grouped"]][["unp_lo", "unp_hi"]].notna().all().all())
check("n_structures >= 1 and <= distinct chains per canonical", (canon["n_structures"] >= 1).all())
check("rigidity vocabulary", set(canon["rigidity"]) <= {"rigid", "flexible", "intermediate", "rigid_helical",
                                                         "flexible_disordered", "strand_spanning", "undetermined"})
_hc = canon[canon["rigidity_confidence"] == "high"]
check("every high-confidence rigidity has >=3 structures or is helical",
      ((_hc["n_structures"] >= 3) | (_hc["rigidity"] == "rigid_helical")).all())
check("dom_dist_cv present whenever rigidity is high-confidence and not purely helical",
      canon[(canon["rigidity_confidence"] == "high") & (canon["rigidity"] != "rigid_helical")]["dom_dist_cv"].notna().all())
check("chimeras excluded from canonical layer",
      not set(canon["uniprot"].dropna()) & {"P00720"} if "P00720" in set(inst["uniprot"].dropna()) else True)
check("length_min <= length_median <= length_max",
      ((canon["length_min"] <= canon["length_median"]) & (canon["length_median"] <= canon["length_max"])).all())
check("liability columns present", all(c in canon for c in ("n_glyc", "deamid_ng", "isomer_dg", "n_cys", "furin")))
if "rigidity_final" in canon:
    check("rigidity_final vocabulary",
          set(canon["rigidity_final"]) <= {"rigid", "flexible", "intermediate", "rigid_helical",
                                           "flexible_disordered", "strand_spanning", "undetermined"})
    check("superposition used only where it could be computed",
          (canon.loc[canon["rigidity_source"] == "superposition", "superpose_rigidity"].notna()
           & (canon.loc[canon["rigidity_source"] == "superposition", "superpose_rigidity"] != "undetermined")).all())
    check("superpose only on >=3 structures",
          (canon.loc[canon["superpose_rigidity"].notna(), "n_structures"] >= 3).all())

print()
print("instances:", len(inst), "qc_pass:", int(ok["qc_pass"].sum()),
      "from", inst.groupby(["pdb_id", "auth_asym_id"]).ngroups, "chains")
print("canonical linkers:", len(canon), " grouped:", int(canon["grouped"].sum()),
      " distinct UniProt:", canon["uniprot"].nunique(), " human:", int(canon["is_human"].sum()))
print("rigidity:", dict(canon["rigidity"].value_counts()))
print("confidence:", dict(canon["rigidity_confidence"].value_counts()))
print(">=3 structures:", int((canon["n_structures"] >= 3).sum()),
      " max structures:", int(canon["n_structures"].max()))
print("carry a hard liability (glyc/NG/DG/furin):",
      int(((canon[["n_glyc", "deamid_ng", "isomer_dg", "furin"]] > 0).any(axis=1)).sum()))
print("run_info:", {k: info[k] for k in ("n_entries", "n_instances", "n_canonical", "n_canonical_grouped") if k in info})
sys.exit(1 if fails else 0)
