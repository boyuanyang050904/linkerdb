"""Validate the superposition-based rigidity measure on the golden set.

For each golden protein, take its best-supported canonical linker, re-read the
structures of that linker's instances, superpose the N-domain and measure how far
the C-domain moves, and compare the superposition call to the expected rigidity and
to the cheap dom_dist/angle-range call.

Usage: python scripts/superpose_golden.py [DB_DIR]
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from linkerdb.config import Settings                      # noqa: E402
from linkerdb.domains.cath import CathProvider, cath_paths  # noqa: E402
from linkerdb.download import Downloader, fetch_structure  # noqa: E402
from linkerdb.golden import coarse, load_golden            # noqa: E402
from linkerdb.linkers import map_domains                   # noqa: E402
from linkerdb.sifts import load_sifts, sifts_path          # noqa: E402
from linkerdb.structure import load_structure              # noqa: E402
from linkerdb.superpose import classify, domain_unp_ca, mobility  # noqa: E402

db = Path(sys.argv[1] if len(sys.argv) > 1 else "results/db")
settings = Settings()
dl = Downloader(timeout=settings.timeout, retries=settings.retries)
p = cath_paths(settings.raw_dir)
provider = CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])
sifts = load_sifts(sifts_path(settings.raw_dir))
canonical = pd.read_parquet(db / "linkers_canonical.parquet")
instances = pd.read_parquet(db / "linkers_instances.parquet")
golden = load_golden()

rows = []
for e in golden.itertuples(index=False):
    m = canonical[canonical["uniprot"] == e.uniprot]
    if m.empty:
        continue
    top = m.sort_values("n_structures", ascending=False).iloc[0]
    sub = instances[(instances["uniprot"] == e.uniprot) & instances["qc_pass"]
                    & (instances["unp_lo"] <= top["unp_hi"]) & (instances["unp_hi"] >= top["unp_lo"])]
    sub = sub.sort_values("resolution")                    # best resolution first = reference
    doms = []
    for r in sub.itertuples():
        try:
            path, fmt = fetch_structure(r.pdb_id, settings, dl)
            model = load_structure(path, r.pdb_id, fmt)
            if r.auth_asym_id not in model.chains:
                continue
            ch = model.chains[r.auth_asym_id]
            ch.use_sifts(sifts.get(r.pdb_id, {}).get(r.auth_asym_id))
            cd = provider.get(r.pdb_id, r.cath_chain)
            mapped = map_domains(ch, cd, settings.rules)
            by_id = {d.domain_id: d for d in mapped.domains}
            if r.domain_n not in by_id or r.domain_c not in by_id:
                continue
            nca = domain_unp_ca(ch, by_id[r.domain_n], e.uniprot)
            cca = domain_unp_ca(ch, by_id[r.domain_c], e.uniprot)
            if nca and cca:
                doms.append((r.pdb_id, nca, cca))
        except Exception as exc:  # noqa: BLE001 - best effort per structure
            print(f"  skip {r.pdb_id}: {exc}", file=sys.stderr)
    mob = mobility(doms)
    rows.append({
        "protein": e.protein, "expected": e.expected_rigidity,
        "proxy": coarse(top["rigidity"]), "proxy_label": top["rigidity"],
        "superpose": classify(mob["max_c_rmsd"], mob["max_hinge"]),
        "max_c_rmsd": None if mob["max_c_rmsd"] is None else round(mob["max_c_rmsd"], 1),
        "max_hinge": None if mob["max_hinge"] is None else round(mob["max_hinge"], 0),
        "n_used": mob["n_used"],
    })

out = pd.DataFrame(rows)
out.to_csv(db / "superpose_golden.csv", index=False)
print(out.to_string(index=False))


def agree(col):
    d = out[~out[col].isin(["intermediate", "undetermined"])]
    return f"{(d[col] == d['expected']).sum()}/{len(d)}"


print(f"\nproxy (dom_dist+angle) agrees: {agree('proxy')}")
print(f"superposition agrees:          {agree('superpose')}")
print(f"wrote {db / 'superpose_golden.csv'}")
