"""Compute superposition-based rigidity for every high-confidence canonical linker
and fold it into results/db/linkers_canonical.parquet.

For each canonical linker with >= MIN_STRUCT structures, re-read its instances'
structures (cached), superpose the N-domain and measure the C-domain displacement and
hinge angle (superpose.py). Adds columns superpose_c_rmsd / superpose_hinge /
superpose_n_used / superpose_rigidity, and sets rigidity_final (superposition where it
could be computed, else the cheap proxy) + rigidity_source.

Structures are grouped so each chain is read once. Usage: python scripts/superpose_db.py [DB_DIR]
"""

import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from linkerdb.aggregate import group_instances                 # noqa: E402
from linkerdb.config import Settings                            # noqa: E402
from linkerdb.domains.cath import CathProvider, cath_paths       # noqa: E402
from linkerdb.download import Downloader, fetch_structure        # noqa: E402
from linkerdb.linkers import map_domains                         # noqa: E402
from linkerdb.sifts import load_sifts, sifts_path                # noqa: E402
from linkerdb.structure import load_structure                    # noqa: E402
from linkerdb.superpose import classify, domain_unp_ca, mobility  # noqa: E402

MIN_STRUCT = 3
db = Path(sys.argv[1] if len(sys.argv) > 1 else "results/db")
settings = Settings()
dl = Downloader(timeout=settings.timeout, retries=settings.retries)
p = cath_paths(settings.raw_dir)
provider = CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])
sifts = load_sifts(sifts_path(settings.raw_dir))
canonical = pd.read_parquet(db / "linkers_canonical.parquet")
instances = pd.read_parquet(db / "linkers_instances.parquet")

# which canonical linkers to do, and which instances belong to each (same grouping as aggregate)
target = set(canonical.loc[canonical["n_structures"] >= MIN_STRUCT, "canonical_id"])
print(f"targets: {len(target)} canonical linkers with >= {MIN_STRUCT} structures", file=sys.stderr)
jobs_by_chain: dict[tuple, list] = defaultdict(list)   # (pdb, auth_chain) -> [(canonical_id, cath_chain, dn, dc, res, uniprot)]
for cid, sub, grouped in group_instances(instances):
    if cid not in target:
        continue
    acc = sub["uniprot"].iloc[0]
    for r in sub.itertuples():
        jobs_by_chain[(r.pdb_id, r.auth_asym_id)].append((cid, r.cath_chain, r.domain_n, r.domain_c,
                                                          r.resolution, acc))

ca_by_cid: dict[str, list] = defaultdict(list)         # cid -> [(resolution, n_ca, c_ca)]
t0 = time.time()
chains = sorted(jobs_by_chain)
for k, (pdb_id, auth) in enumerate(chains):
    if k % 2000 == 0:
        print(f"  {k}/{len(chains)} chains, {time.time() - t0:.0f}s", file=sys.stderr)
    try:
        path, fmt = fetch_structure(pdb_id, settings, dl)
        model = load_structure(path, pdb_id, fmt)
        if auth not in model.chains:
            continue
        ch = model.chains[auth]
        ch.use_sifts(sifts.get(pdb_id, {}).get(auth))
        mapped_cache = {}
        for cid, cath_chain, dn, dc, res, acc in jobs_by_chain[(pdb_id, auth)]:
            if cath_chain not in mapped_cache:
                try:
                    cd = provider.get(pdb_id, cath_chain)
                    mapped_cache[cath_chain] = {d.domain_id: d for d in map_domains(ch, cd, settings.rules).domains}
                except Exception:  # noqa: BLE001
                    mapped_cache[cath_chain] = {}
            by_id = mapped_cache[cath_chain]
            if dn in by_id and dc in by_id:
                nca = domain_unp_ca(ch, by_id[dn], acc)
                cca = domain_unp_ca(ch, by_id[dc], acc)
                if nca and cca:
                    ca_by_cid[cid].append((res if pd.notna(res) else 99.0, nca, cca))
    except Exception as exc:  # noqa: BLE001
        print(f"  skip {pdb_id}_{auth}: {exc}", file=sys.stderr)

# compute mobility per canonical linker
rec = {}
for cid, items in ca_by_cid.items():
    items.sort(key=lambda x: x[0])                     # best resolution first = reference
    mob = mobility([(f"s{i}", n, c) for i, (_, n, c) in enumerate(items)])
    rec[cid] = (mob["max_c_rmsd"], mob["max_hinge"], mob["n_used"],
                classify(mob["max_c_rmsd"], mob["max_hinge"]))

canonical["superpose_c_rmsd"] = canonical["canonical_id"].map(lambda c: rec.get(c, (None,) * 4)[0])
canonical["superpose_hinge"] = canonical["canonical_id"].map(lambda c: rec.get(c, (None,) * 4)[1])
canonical["superpose_n_used"] = canonical["canonical_id"].map(lambda c: rec.get(c, (None,) * 4)[2])
canonical["superpose_rigidity"] = canonical["canonical_id"].map(lambda c: rec.get(c, (None,) * 4)[3])

usable = canonical["superpose_rigidity"].notna() & (canonical["superpose_rigidity"] != "undetermined")
canonical["rigidity_final"] = canonical["rigidity"].where(~usable, canonical["superpose_rigidity"])
canonical["rigidity_source"] = usable.map({True: "superposition", False: "proxy"})
canonical.to_parquet(db / "linkers_canonical.parquet", index=False)

print(f"\ndone in {time.time() - t0:.0f}s; superposition computed for {int(usable.sum())} linkers "
      f"({canonical['superpose_rigidity'].eq('undetermined').sum()} undetermined, rest <3 structures)")
print("rigidity_final:", dict(canonical["rigidity_final"].value_counts()))
print("changed vs proxy:", int((canonical["rigidity_final"] != canonical["rigidity"]).sum()))
