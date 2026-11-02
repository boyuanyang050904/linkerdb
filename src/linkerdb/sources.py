"""Entry-level metadata: resolution, experiment type, entity mapping, sequence clusters."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterable

from .config import ENTRY_TYPE_URL, RCSB_CLUSTERS_URL, RCSB_GRAPHQL_URL, RESOLU_URL
from .download import Downloader, DownloadError

log = logging.getLogger(__name__)


def index_paths(raw_dir: Path, identity: int = 40) -> dict[str, Path]:
    d = Path(raw_dir) / "index"
    return {"resolu": d / "resolu.idx", "entry_type": d / "pdb_entry_type.txt",
            "clusters": d / f"clusters-by-entity-{identity}.txt"}


def fetch_indices(raw_dir: Path, downloader: Downloader, identity: int = 40) -> dict[str, Path]:
    paths = index_paths(raw_dir, identity)
    downloader.fetch(RESOLU_URL, paths["resolu"])
    downloader.fetch(ENTRY_TYPE_URL, paths["entry_type"])
    downloader.fetch(RCSB_CLUSTERS_URL.format(pct=identity), paths["clusters"])
    return paths


def parse_resolu(path: Path) -> dict[str, float | None]:
    """``resolu.idx`` -> ``{pdb_id: resolution or None}`` (``-1.00`` = no resolution, e.g. NMR)."""
    out: dict[str, float | None] = {}
    with open(path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split(";")]
            if len(parts) != 2 or len(parts[0]) != 4:
                continue
            try:
                val = float(parts[1])
            except ValueError:
                continue
            out[parts[0].lower()] = val if val > 0 else None
    return out


def parse_entry_type(path: Path) -> dict[str, tuple[str, str]]:
    """``pdb_entry_type.txt`` -> ``{pdb_id: (molecule type, method class)}``."""
    out = {}
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) >= 3:
                out[t[0].lower()] = (t[1], t[2])
    return out


def parse_clusters(path: Path, wanted: set[str] | None = None) -> dict[str, int]:
    """``clusters-by-entity-NN.txt`` -> ``{"1ABC_1": cluster_index}``."""
    out = {}
    with open(path) as fh:
        for k, line in enumerate(fh):
            for ent in line.split():
                if wanted is None or ent in wanted:
                    out[ent] = k
    return out


GRAPHQL_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_info { resolution_combined deposited_model_count }
    exptl { method }
    polymer_entities {
      rcsb_id
      entity_poly { type }
      rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids uniprot_ids }
    }
  }
}
"""


def _flatten_entry(e: dict) -> dict:
    ents = []
    for pe in e.get("polymer_entities") or []:
        ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        ents.append({
            "entity": pe["rcsb_id"],
            "type": (pe.get("entity_poly") or {}).get("type"),
            "auth_asym_ids": ids.get("auth_asym_ids") or [],
            "uniprot_ids": ids.get("uniprot_ids") or [],
        })
    info = e.get("rcsb_entry_info") or {}
    res = info.get("resolution_combined") or []
    return {
        "pdb_id": e["rcsb_id"].lower(),
        "methods": [x["method"] for x in e.get("exptl") or []],
        "resolution": min(res) if res else None,
        "n_models": info.get("deposited_model_count"),
        "entities": ents,
    }


def fetch_entity_info(pdb_ids: Iterable[str], cache: Path, downloader: Downloader,
                      batch: int = 400) -> dict[str, dict]:
    """Entity/chain/UniProt metadata from RCSB GraphQL, cached as JSON lines.

    Resumable: entries already present in ``cache`` are not queried again.
    """
    cache = Path(cache)
    cache.parent.mkdir(parents=True, exist_ok=True)
    have: dict[str, dict] = {}
    if cache.exists():
        with open(cache) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue  # truncated last line from an interrupted run
                have[rec["pdb_id"]] = rec
    todo = sorted({p.lower() for p in pdb_ids} - set(have))
    with open(cache, "a") as out:
        for k in range(0, len(todo), batch):
            ids = [p.upper() for p in todo[k:k + batch]]
            try:
                data = downloader.get_json(RCSB_GRAPHQL_URL, method="POST",
                                           payload={"query": GRAPHQL_QUERY, "variables": {"ids": ids}})
            except DownloadError as exc:
                log.warning("GraphQL batch %d failed: %s", k // batch, exc)
                continue
            got = set()
            for e in (data.get("data") or {}).get("entries") or []:
                if e is None:
                    continue
                rec = _flatten_entry(e)
                have[rec["pdb_id"]] = rec
                got.add(rec["pdb_id"])
                out.write(json.dumps(rec) + "\n")
            for missing in {i.lower() for i in ids} - got:  # obsolete / withdrawn entries
                rec = {"pdb_id": missing, "missing": True, "entities": []}
                have[missing] = rec
                out.write(json.dumps(rec) + "\n")
            out.flush()
            log.info("GraphQL: %d / %d entries", min(k + batch, len(todo)), len(todo))
    return have
