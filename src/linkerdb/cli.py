"""Command-line interface.

    linkerdb fetch-meta                 # CATH files + wwPDB/RCSB indices
    linkerdb sample  --out results/sample50
    linkerdb build   --out results/full --max-resolution 2.5
    linkerdb summary results/full/linkers.parquet
    linkerdb show 13pk_A_1
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import random
import sys
import time
from pathlib import Path

import pandas as pd

from . import __version__
from .aggregate import aggregate
from .config import Rules, Settings
from .domains.resolve import SOURCES, build_provider
from .download import Downloader, fetch_structure
from .linkers import ChainSkip, chain_ss, extract_linkers, map_domains, natural_problem
from .pipeline import Run, _json_default, write_outputs
from .retrieve import format_hits, query_canonical
from .selection import (attach_clusters, candidate_chains, chains_of_entries, dedup_linkers, rank_candidates,
                        sample_entries)
from .sifts import fetch_sifts, load_sifts
from .sources import fetch_entity_info, fetch_indices, parse_clusters, parse_entry_type, parse_resolu
from .stats import show_linker, summarize, summarize_canonical
from .structure import load_structure

log = logging.getLogger("linkerdb")


def _settings(args) -> Settings:
    rules = Rules(**{k: v for k, v in {
        "max_ext": getattr(args, "max_ext", None),
        "max_resolution": getattr(args, "max_resolution", None),
    }.items() if v is not None})
    kw = {"rules": rules, "workers": args.workers}
    if args.data_dir:
        kw["data_dir"] = Path(args.data_dir)
    return Settings(**kw)


def _provider(settings: Settings, downloader: Downloader, order: str):
    prov, notes = build_provider(settings.raw_dir, downloader, tuple(order.split(",")))
    for n in notes:
        log.info("domain source: %s", n)
    return prov, notes


def _jobs(provider, entry_chains: dict[str, list[str]]):
    jobs = []
    for pdb, chains in sorted(entry_chains.items()):
        items = [(c, provider.get(pdb, c)) for c in sorted(chains)]
        items = [(c, cd) for c, cd in items if cd is not None]
        if items:
            jobs.append((pdb, items))
    return jobs


def _info(settings: Settings, provider, notes, **extra) -> dict:
    return {"linkerdb_version": __version__, "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
            "domain_source": provider.source, "domain_source_notes": notes,
            "rules": dataclasses.asdict(settings.rules), **extra}


def cmd_fetch_meta(args) -> int:
    settings = _settings(args)
    dl = Downloader(timeout=settings.timeout, retries=settings.retries)
    _provider(settings, dl, args.sources)
    paths = fetch_indices(settings.raw_dir, dl)
    for k, p in paths.items():
        print(f"{k}: {p} ({p.stat().st_size / 1e6:.1f} MB)")
    return 0


def cmd_sample(args) -> int:
    settings = _settings(args)
    dl = Downloader(timeout=settings.timeout, retries=settings.retries)
    provider, notes = _provider(settings, dl, args.sources)
    paths = fetch_indices(settings.raw_dir, dl)
    etype = parse_entry_type(paths["entry_type"])
    entries = sample_entries(provider, etype, args.n_diffraction, args.n_nmr, args.seed)
    jobs = _jobs(provider, chains_of_entries(provider, entries))
    sifts = load_sifts(fetch_sifts(settings.raw_dir, dl), set(entries))
    run = Run(settings, args.run_name)
    chains, linkers = run.execute(jobs, workers=args.workers, one_per_entity=True, force=args.force, sifts=sifts)
    info = _info(settings, provider, notes, entries=entries, seed=args.seed)
    final = linkers[linkers["qc_pass"]].reset_index(drop=True) if not linkers.empty else linkers
    out = write_outputs(Path(args.out), chains, linkers, final, info)
    bg = chains.loc[chains["status"] == "ok", "sequence"] if "sequence" in chains else None
    print(summarize(chains, linkers, bg, title=f"sample: {len(entries)} entries"))
    for k, p in out.items():
        print(f"wrote {k}: {p}")
    return 0


def cmd_build(args) -> int:
    """Build the linker database: every domain junction in every multi-domain chain
    (one chain per entity per entry), grouped into canonical linkers by UniProt region."""
    settings = _settings(args)
    rules = settings.rules
    dl = Downloader(timeout=settings.timeout, retries=settings.retries)
    provider, notes = _provider(settings, dl, args.sources)
    paths = fetch_indices(settings.raw_dir, dl, 40)
    resolu = parse_resolu(paths["resolu"])
    etype = parse_entry_type(paths["entry_type"])
    cands = candidate_chains(provider, resolu, etype, rules.max_resolution)
    log.info("candidates: %d chains in %d entries (resolution <= %.2f)", len(cands), cands["pdb_id"].nunique(),
             rules.max_resolution)
    ent_info = fetch_entity_info(cands["pdb_id"].unique(), settings.raw_dir / "rcsb" / "entities.jsonl", dl)
    wanted = {e["entity"] for p in cands["pdb_id"].unique() for e in ent_info.get(p, {}).get("entities", [])}
    clusters = parse_clusters(paths["clusters"], wanted)
    n_dom = {k: provider.get(k[:4], k[4:]).n_domains for k in (cands["pdb_id"] + cands["chain"])}
    cands = attach_clusters(cands, ent_info, clusters, n_dom)
    cands = cands[cands["polymer_type"].isna() | (cands["polymer_type"] == "polypeptide(L)")].reset_index(drop=True)

    if args.sample_clusters:
        nstruct = cands.groupby("cluster")["pdb_id"].nunique()
        pool = sorted(nstruct[(nstruct >= 3) & (nstruct <= 15)].index)  # redundant but cheap to download
        random.Random(args.seed).shuffle(pool)
        pick = set(pool[: args.sample_clusters])
        cands = cands[cands["cluster"].isin(pick)]
        log.info("sample-clusters: %d redundant clusters -> %d chains", len(pick), len(cands))
    elif args.limit:
        keep = sorted(cands["pdb_id"].unique())[: args.limit]
        cands = cands[cands["pdb_id"].isin(keep)]

    cluster_of = {(r.pdb_id, r.chain): r.cluster for r in cands.itertuples()}
    sifts = load_sifts(fetch_sifts(settings.raw_dir, dl), set(cands["pdb_id"]))
    run = Run(settings, args.run_name)
    jobs = _jobs(provider, cands.groupby("pdb_id")["chain"].apply(list).to_dict())
    # one chain per entity per entry: crystallographic copies are not independent evidence
    chains, instances = run.execute(jobs, workers=args.workers, one_per_entity=True, sifts=sifts, force=args.force)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cands.to_parquet(out_dir / "candidates.parquet", index=False)
    if not instances.empty:
        auth_to_cath = {(r.pdb_id, r.auth_asym_id): r.chain for r in chains.itertuples()
                        if isinstance(r.auth_asym_id, str)}
        instances["cath_chain"] = [auth_to_cath.get((p, a), a)
                                   for p, a in zip(instances["pdb_id"], instances["auth_asym_id"])]
        instances["cluster40"] = [cluster_of.get((p, c)) for p, c in zip(instances["pdb_id"], instances["cath_chain"])]
        canonical = aggregate(instances)
    else:
        canonical = pd.DataFrame()
    info = _info(settings, provider, notes, mode="database", n_candidate_chains=len(cands),
                 n_entries=int(cands["pdb_id"].nunique()),
                 n_instances=0 if instances.empty else int(instances["qc_pass"].sum()),
                 n_canonical=len(canonical),
                 n_canonical_grouped=0 if canonical.empty else int(canonical["grouped"].sum()))
    chains.to_parquet(out_dir / "chains.parquet", index=False)
    instances.to_parquet(out_dir / "linkers_instances.parquet", index=False)
    canonical.to_parquet(out_dir / "linkers_canonical.parquet", index=False)
    (out_dir / "run_info.json").write_text(json.dumps(info, indent=2, default=_json_default))
    print(summarize_canonical(canonical, instances, title=f"linker database (<= {rules.max_resolution} Å)"))
    for name in ("chains", "linkers_instances", "linkers_canonical"):
        print(f"wrote {name}: {out_dir / (name + '.parquet')}")
    return 0


def cmd_query(args) -> int:
    df = pd.read_parquet(args.canonical)
    hits = query_canonical(df, distance=args.distance, distance_tol=args.distance_tol, rigidity=args.rigidity,
                           human_only=args.human, min_length=args.min_length, max_length=args.max_length,
                           min_structures=args.min_structures, no_glyc=args.no_glyc, no_cys=args.no_cys,
                           max_liability=args.max_liability, top=args.top)
    print(format_hits(hits))
    return 0


def cmd_summary(args) -> int:
    linkers = pd.read_parquet(args.linkers)
    chains_path = Path(args.linkers).with_name("chains.parquet")
    chains = pd.read_parquet(chains_path) if chains_path.exists() else pd.DataFrame()
    print(summarize(chains, linkers, title=str(args.linkers)))
    return 0


def cmd_show(args) -> int:
    settings = _settings(args)
    dl = Downloader(timeout=settings.timeout, retries=settings.retries)
    provider, _ = _provider(settings, dl, args.sources)
    for lid in args.linker_ids:
        pdb, chain, idx = lid.rsplit("_", 2)
        path, fmt = fetch_structure(pdb, settings, dl)
        model = load_structure(path, pdb, fmt)
        cd = provider.get(pdb, chain)
        if cd is None or chain not in model.chains:
            print(f"{lid}: chain not in domain source / structure")
            continue
        m = model.chains[chain]
        m.use_sifts(load_sifts(fetch_sifts(settings.raw_dir, dl), {pdb}).get(pdb, {}).get(chain))
        try:
            mapped = map_domains(m, cd, settings.rules)
        except ChainSkip as skip:
            print(f"{lid}: {skip}")
            continue
        ss = chain_ss(m)
        nat = natural_problem(m)
        recs = extract_linkers(m, mapped, ss, settings.rules, [nat[0]] if nat else [])
        rec = next((r for r in recs if r["junction_index"] == int(idx)), None)
        if rec is None:
            print(f"{lid}: no such junction (chain has {len(recs)})")
            continue
        print(f"{m.description} | {m.entry.method} {m.entry.resolution} | UniProt {','.join(m.uniprot_accessions())}")
        print(show_linker(m, rec, ss, args.pad))
        if args.json:
            print(json.dumps({k: v for k, v in rec.items() if k not in ("residues_3",)}, default=str, indent=1))
        print()
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="linkerdb", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", help="cache directory (default: $LINKERDB_DATA or ./data)")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--sources", default=",".join(SOURCES), help="domain source preference (default: cath,scope,pfam)")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch-meta", help="download domain and index files")
    s.set_defaults(func=cmd_fetch_meta)

    s = sub.add_parser("sample", help="run the pipeline on a random sample")
    s.add_argument("--n-diffraction", type=int, default=40)
    s.add_argument("--n-nmr", type=int, default=10)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--max-ext", type=int)
    s.add_argument("--run-name", default="sample50")
    s.add_argument("--out", default="results/sample50")
    s.add_argument("--force", action="store_true", help="ignore checkpoints")
    s.set_defaults(func=cmd_sample)

    s = sub.add_parser("build", help="build the linker database (all multi-domain chains -> canonical linkers)")
    s.add_argument("--max-resolution", type=float)
    s.add_argument("--max-ext", type=int)
    s.add_argument("--limit", type=int, help="only the first N candidate entries (testing)")
    s.add_argument("--sample-clusters", type=int,
                   help="process all members of N random redundant clusters (fast, redundancy-rich testing)")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--run-name", default="db")
    s.add_argument("--out", default="results/db")
    s.add_argument("--force", action="store_true", help="ignore checkpoints")
    s.set_defaults(func=cmd_build)

    s = sub.add_parser("query", help="search canonical linkers (e.g. --distance 25 --rigidity rigid --human)")
    s.add_argument("--canonical", default="results/db/linkers_canonical.parquet")
    s.add_argument("--distance", type=float, help="target inter-domain COM distance (Å)")
    s.add_argument("--distance-tol", type=float, default=5.0)
    s.add_argument("--rigidity", choices=["rigid", "flexible", "any"])
    s.add_argument("--human", action="store_true")
    s.add_argument("--min-length", type=int)
    s.add_argument("--max-length", type=int)
    s.add_argument("--min-structures", type=int, default=1)
    s.add_argument("--no-glyc", action="store_true", help="exclude linkers with an N-glycosylation sequon")
    s.add_argument("--no-cys", action="store_true")
    s.add_argument("--max-liability", type=int)
    s.add_argument("--top", type=int, default=20)
    s.set_defaults(func=cmd_query)

    s = sub.add_parser("summary", help="print statistics of a linkers parquet file")
    s.add_argument("linkers")
    s.set_defaults(func=cmd_summary)

    s = sub.add_parser("show", help="print a readable view of linkers, e.g. 13pk_A_1")
    s.add_argument("linker_ids", nargs="+")
    s.add_argument("--pad", type=int, default=15)
    s.add_argument("--json", action="store_true")
    s.set_defaults(func=cmd_show)

    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s", stream=sys.stderr)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
