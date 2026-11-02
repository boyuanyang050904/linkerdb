"""Orchestration: download -> parse -> map domains -> extract linkers, with
per-chain checkpoints so that interrupted runs resume where they stopped."""

from __future__ import annotations

import json
import logging
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .config import Rules, Settings
from .domains.base import ChainDomains
from .download import Downloader, NotFoundError, fetch_structures
from .linkers import ChainSkip, chain_ss, extract_linkers, map_domains, natural_problem
from .structure import load_structure

log = logging.getLogger(__name__)

Job = tuple[str, list[tuple[str, ChainDomains]]]  # (pdb_id, [(cath chain, domains), ...])


# ---------------------------------------------------------------- one entry
def _chain_record(pdb_id: str, chain: str, cd: ChainDomains, status: str, detail: str = "", **kw) -> dict:
    rec = {"pdb_id": pdb_id, "chain": chain, "domain_source": cd.source, "n_domains": cd.n_domains,
           "status": status, "detail": detail}
    rec.update(kw)
    return rec


def _ss_agreement(ss: str, annotated: str) -> float | None:
    pairs = [(a, b) for a, b in zip(ss, annotated) if a != "?"]
    if not pairs or set(annotated) == {"-"}:
        return None
    three = lambda c: "H" if c in "HGI" else ("E" if c == "E" else "-")  # noqa: E731
    return float(np.mean([three(a) == three(b) for a, b in pairs]))


def process_entry(pdb_id: str, chains: Sequence[tuple[str, ChainDomains]], path: Path | None, fmt: str,
                  rules: Rules, one_per_entity: bool = False, fetch_error: str | None = None,
                  unp: dict | None = None) -> dict:
    """Process the requested chains of one entry; never raises.

    ``unp`` optionally maps author chain -> SIFTS UniProt segments (preferred over
    the mmCIF ``_struct_ref_seq`` mapping).
    """
    out = {"pdb_id": pdb_id, "chains": [], "linkers": []}
    if path is None:
        for ch, cd in chains:
            out["chains"].append(_chain_record(pdb_id, ch, cd, "download_failed", fetch_error or ""))
        return out
    try:
        model = load_structure(path, pdb_id, fmt)
    except Exception as exc:  # noqa: BLE001 - any parser failure is recorded
        for ch, cd in chains:
            out["chains"].append(_chain_record(pdb_id, ch, cd, "parse_failed", f"{type(exc).__name__}: {exc}"))
        return out

    resolved = []
    for ch, cd in chains:
        auth, extra = ch, []
        if auth not in model.chains:
            prot = [c for c, m in model.chains.items() if m.polymer_type == "polypeptide(L)"]
            if ch == "0" and len(prot) == 1:
                auth, extra = prot[0], ["chain0_remap"]
            else:
                out["chains"].append(_chain_record(pdb_id, ch, cd, "chain_not_found"))
                continue
        resolved.append((ch, auth, cd, extra))

    if one_per_entity:
        # identical copies of one entity (homo-oligomers) give identical linkers: keep the most complete
        best: dict[str, tuple] = {}
        for item in resolved:
            m = model.chains[item[1]]
            score = (int(m.observed.sum()), -ord(item[1][0]))
            if m.entity_id not in best or score > best[m.entity_id][0]:
                best[m.entity_id] = (score, item)
        keep = {id(v[1]) for v in best.values()}
        for item in resolved:
            if id(item) not in keep:
                out["chains"].append(_chain_record(pdb_id, item[0], item[2], "duplicate_entity_copy",
                                                   entity_id=model.chains[item[1]].entity_id))
        resolved = [item for item in resolved if id(item) in keep]

    for ch, auth, cd, extra in resolved:
        m = model.chains[auth]
        if unp:
            m.use_sifts(unp.get(auth))
        base = dict(auth_asym_id=auth, entity_id=m.entity_id, method=m.entry.method,
                    resolution=m.entry.resolution, n_models=m.entry.n_models, chain_length=m.n,
                    n_observed=int(m.observed.sum()), uniprot=",".join(m.uniprot_accessions()),
                    description=m.description, taxid=m.taxid, fmt=m.entry.fmt, sequence=m.sequence,
                    unp_source=m.ref_source)
        nat = natural_problem(m)
        flags = list(extra)
        if nat:
            flags.append(nat[0])
        try:
            mapped = map_domains(m, cd, rules)
        except ChainSkip as skip:
            out["chains"].append(_chain_record(pdb_id, ch, cd, skip.reason, skip.detail, **base))
            continue
        ss = chain_ss(m)
        linkers = extract_linkers(m, mapped, ss, rules, extra_flags=flags)
        for rec in linkers:
            rec["domain_source"] = cd.source
            rec["description"] = m.description
            rec["taxid"] = m.taxid
        status = f"excluded_{nat[0]}" if nat else "ok"
        out["chains"].append(_chain_record(
            pdb_id, ch, cd, status, nat[1] if nat else "", numbering_check=mapped.numbering_check,
            chain_flags=sorted(set(mapped.flags + flags)), n_junctions=len(linkers),
            n_linkers_pass=sum(r["qc_pass"] for r in linkers), dssp_vs_annotation=_ss_agreement(ss, m.annotated_ss),
            **base))
        out["linkers"].extend(linkers)
    return out


def _worker(args) -> dict:
    pdb_id, chains, path, fmt, rules, one_per_entity, err, unp = args
    return process_entry(pdb_id, chains, path, fmt, rules, one_per_entity, err, unp)


# ---------------------------------------------------------------- runs
class Run:
    """A named run directory with per-chain checkpoints."""

    def __init__(self, settings: Settings, name: str) -> None:
        self.settings = settings
        self.name = name
        self.dir = settings.work_dir / name
        self.ckpt = self.dir / "entries"
        self.ckpt.mkdir(parents=True, exist_ok=True)

    def _ckpt_path(self, pdb_id: str, chain: str) -> Path:
        return self.ckpt / f"{pdb_id.lower()}_{chain}.json"

    def done(self, pdb_id: str, chain: str) -> bool:
        return self._ckpt_path(pdb_id, chain).exists()

    def load(self, pdb_id: str, chain: str) -> dict:
        return json.loads(self._ckpt_path(pdb_id, chain).read_text())

    def save(self, result: dict) -> None:
        by_chain: dict[str, dict] = {}
        for c in result["chains"]:
            by_chain[c["chain"]] = {"chains": [c], "linkers": []}
        for rec in result["linkers"]:
            # linker records carry the author chain; map back to the requested (CATH) chain
            for c in result["chains"]:
                if c.get("auth_asym_id") == rec["auth_asym_id"]:
                    by_chain[c["chain"]]["linkers"].append(rec)
                    break
        for ch, payload in by_chain.items():
            tmp = self._ckpt_path(result["pdb_id"], ch).with_suffix(".tmp")
            tmp.write_text(json.dumps(payload, default=_json_default))
            tmp.replace(self._ckpt_path(result["pdb_id"], ch))

    def execute(self, jobs: Iterable[Job], workers: int | None = None, one_per_entity: bool = False,
                downloader: Downloader | None = None, force: bool = False,
                sifts: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Run (or resume) ``jobs`` and return ``(chains_df, linkers_df)`` for exactly those jobs."""
        jobs = [(p.lower(), list(chs)) for p, chs in jobs]
        todo = []
        for pdb_id, chs in jobs:
            pending = [(c, cd) for c, cd in chs if force or not self.done(pdb_id, c)]
            if pending:
                # entity de-duplication needs all chains of the entry together
                todo.append((pdb_id, chs if one_per_entity else pending))
        log.info("%s: %d entries total, %d to process", self.name, len(jobs), len(todo))
        if todo:
            t0 = time.time()
            fetched = fetch_structures([p for p, _ in todo], self.settings, downloader)
            log.info("downloaded/cached %d entries in %.1fs", len(fetched), time.time() - t0)
            args = []
            for pdb_id, chs in todo:
                res = fetched.get(pdb_id.upper())
                if isinstance(res, Exception):
                    status = "obsolete_or_missing" if isinstance(res, NotFoundError) else "download"
                    args.append((pdb_id, chs, None, "cif", self.settings.rules, one_per_entity, f"{status}: {res}",
                                 None))
                else:
                    args.append((pdb_id, chs, res[0], res[1], self.settings.rules, one_per_entity, None,
                                 (sifts or {}).get(pdb_id)))
            n_workers = workers or self.settings.workers
            done = 0
            if n_workers > 1:
                # forkserver: never fork a process that has run download threads
                ctx = multiprocessing.get_context("forkserver")
                with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
                    for result in pool.map(_worker, args, chunksize=4):
                        self.save(result)
                        done += 1
                        if done % 200 == 0:
                            log.info("processed %d / %d entries", done, len(args))
            else:
                for a in args:
                    self.save(_worker(a))
            log.info("processed %d entries in %.1fs", len(args), time.time() - t0)
        chains, linkers = [], []
        for pdb_id, chs in jobs:
            for c, _ in chs:
                if not self.done(pdb_id, c):
                    continue
                payload = self.load(pdb_id, c)
                chains.extend(payload["chains"])
                linkers.extend(payload["linkers"])
        return pd.DataFrame(chains), pd.DataFrame(linkers)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def write_outputs(out_dir: Path, chains: pd.DataFrame, linkers: pd.DataFrame,
                  final: pd.DataFrame | None, info: dict) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {"chains": out_dir / "chains.parquet", "linkers_all": out_dir / "linkers_all.parquet"}
    chains.to_parquet(paths["chains"], index=False)
    linkers.to_parquet(paths["linkers_all"], index=False)
    if final is not None:
        paths["linkers"] = out_dir / "linkers.parquet"
        final.to_parquet(paths["linkers"], index=False)
    info = dict(info)
    info["rules"] = info.get("rules") or {}
    (out_dir / "run_info.json").write_text(json.dumps(info, indent=2, default=_json_default))
    return paths


def rules_dict(rules: Rules) -> dict:
    return asdict(rules)
