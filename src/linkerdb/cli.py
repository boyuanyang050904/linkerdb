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
import sys
import time
from pathlib import Path

import pandas as pd

from . import __version__
from .config import Rules, Settings
from .domains.resolve import SOURCES, build_provider
from .download import Downloader, fetch_structure
from .linkers import ChainSkip, chain_ss, extract_linkers, map_domains, natural_problem
from .pipeline import Run, write_outputs
from .selection import (attach_clusters, candidate_chains, chains_of_entries, dedup_linkers, rank_candidates,
                        sample_entries)
from .sifts import fetch_sifts, load_sifts
from .sources import fetch_entity_info, fetch_indices, parse_clusters, parse_entry_type, parse_resolu
from .stats import show_linker, summarize
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
    settings = _settings(args)
    rules = settings.rules
    dl = Downloader(timeout=settings.timeout, retries=settings.retries)
    provider, notes = _provider(settings, dl, args.sources)
    paths = fetch_indices(settings.raw_dir, dl, args.identity)
    resolu = parse_resolu(paths["resolu"])
    etype = parse_entry_type(paths["entry_type"])
    cands = candidate_chains(provider, resolu, etype, rules.max_resolution)
    log.info("candidates: %d chains in %d entries (resolution <= %.2f)", len(cands), cands["pdb_id"].nunique(),
             rules.max_resolution)
    if args.limit:
        keep = sorted(cands["pdb_id"].unique())[: args.limit]
        cands = cands[cands["pdb_id"].isin(keep)]
    ent_info = fetch_entity_info(cands["pdb_id"].unique(), settings.raw_dir / "rcsb" / "entities.jsonl", dl)
    wanted = {e["entity"] for p in cands["pdb_id"].unique() for e in ent_info.get(p, {}).get("entities", [])}
    clusters = parse_clusters(paths["clusters"], wanted)
    n_dom = {k: provider.get(k[:4], k[4:]).n_domains for k in (cands["pdb_id"] + cands["chain"])}
    ranked = rank_candidates(attach_clusters(cands, ent_info, clusters, n_dom))
    log.info("clusters: %d (%d ranked chains)", ranked["cluster"].nunique(), len(ranked))

    sifts = load_sifts(fetch_sifts(settings.raw_dir, dl), set(ranked["pdb_id"]))
    run = Run(settings, args.run_name)
    resolved: dict[str, tuple[str, str]] = {}
    all_chains, all_linkers = [], []
    for attempt in range(rules.tries_per_cluster):
        todo = ranked[(ranked["rank"] == attempt) & ~ranked["cluster"].isin(resolved)]
        if todo.empty:
            break
        log.info("attempt %d: %d clusters", attempt + 1, len(todo))
        jobs = _jobs(provider, todo.groupby("pdb_id")["chain"].apply(list).to_dict())
        chains, linkers = run.execute(jobs, workers=args.workers, sifts=sifts)
        all_chains.append(chains)
        all_linkers.append(linkers)
        ok = chains[chains["status"] == "ok"]
        okset = set(zip(ok["pdb_id"], ok["chain"]))
        for row in todo.itertuples():
            if (row.pdb_id, row.chain) in okset:
                resolved[row.cluster] = (row.pdb_id, row.chain)
    chains = pd.concat(all_chains, ignore_index=True) if all_chains else pd.DataFrame()
    linkers = pd.concat(all_linkers, ignore_index=True) if all_linkers else pd.DataFrame()
    rep = {pc: clu for clu, pc in resolved.items()}
    chains["cluster"] = [rep.get((p, c)) for p, c in zip(chains["pdb_id"], chains["chain"])]
    chains["representative"] = chains["cluster"].notna()
    rank_of = {(r.pdb_id, r.chain): (r.cluster, r.rank) for r in ranked.itertuples()}
    if not linkers.empty:
        # linker records carry the author chain; representatives are keyed by the CATH chain
        auth_to_cath = {(r.pdb_id, r.auth_asym_id): r.chain for r in chains.itertuples()
                        if isinstance(r.auth_asym_id, str)}
        linkers["cath_chain"] = [auth_to_cath.get((p, a), a) for p, a in zip(linkers["pdb_id"], linkers["auth_asym_id"])]
        linkers["cluster"] = [rank_of.get((p, c), (None, None))[0] for p, c in zip(linkers["pdb_id"], linkers["cath_chain"])]
        is_rep = [(p, c) in rep for p, c in zip(linkers["pdb_id"], linkers["cath_chain"])]
        final = linkers[pd.Series(is_rep, index=linkers.index) & linkers["qc_pass"]]
        seqs = {(r.pdb_id, r.auth_asym_id): r.sequence for r in chains.itertuples()
                if isinstance(getattr(r, "sequence", None), str)}
        final = dedup_linkers(final, ["resolution", "pdb_id", "auth_asym_id", "junction_index"],
                              rules.dedup_identity, chain_seqs=seqs)
        dropped = final.attrs.get("dropped", {})
        linkers["dedup_dropped"] = linkers["linker_id"].map(dropped)
    else:
        final, dropped = linkers, {}
    info = _info(settings, provider, notes, identity=args.identity, n_candidates=len(cands),
                 n_clusters=int(ranked["cluster"].nunique()), n_representatives=len(resolved),
                 n_dedup_dropped=len(dropped), n_final=len(final))
    Path(args.out).mkdir(parents=True, exist_ok=True)
    ranked.to_parquet(Path(args.out) / "candidates.parquet", index=False)
    out = write_outputs(Path(args.out), chains, linkers, final, info)
    rep_chains = chains[chains["representative"]]
    bg = rep_chains["sequence"] if "sequence" in rep_chains else None
    print(summarize(rep_chains, final, bg, title=f"full dataset (<= {rules.max_resolution} Å, {args.identity}% clusters)"))
    print(f"clusters: {ranked['cluster'].nunique()}  representatives: {len(resolved)}  "
          f"dedup dropped: {len(dropped)}  final linkers: {len(final)}")
    for k, p in out.items():
        print(f"wrote {k}: {p}")
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

    s = sub.add_parser("build", help="full non-redundant dataset")
    s.add_argument("--max-resolution", type=float)
    s.add_argument("--identity", type=int, default=40, choices=[30, 40, 50, 70, 90, 95, 100])
    s.add_argument("--max-ext", type=int)
    s.add_argument("--limit", type=int, help="only the first N candidate entries (testing)")
    s.add_argument("--run-name", default="full")
    s.add_argument("--out", default="results/full")
    s.set_defaults(func=cmd_build)

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
