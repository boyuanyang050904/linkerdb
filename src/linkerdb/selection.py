"""Dataset selection: the 50-structure sample, full-run candidates and cluster
representatives (docs/design.md §8).

Linker-level redundancy removal (sequence alignment + k-mer pre-filter) was removed
in v3: the database groups instances by UniProt region (see ``aggregate.py``) instead
of deleting redundant linkers, so the pairwise-alignment machinery is no longer needed.
"""

from __future__ import annotations

import logging
import random
from collections import defaultdict

import pandas as pd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------- sample
def sample_entries(provider, entry_type: dict[str, tuple[str, str]], n_diffraction: int = 40, n_nmr: int = 10,
                   seed: int = 0) -> list[str]:
    """Random multi-domain entries: ``n_diffraction`` diffraction (X-ray/neutron) + ``n_nmr`` NMR entries."""
    entries = sorted(provider.entries())
    prot = [e for e in entries if entry_type.get(e, ("", ""))[0] in ("prot", "prot-nuc")]
    diff = [e for e in prot if entry_type[e][1] == "diffraction"]
    nmr = [e for e in prot if entry_type[e][1] == "NMR"]
    rng = random.Random(seed)
    return sorted(rng.sample(diff, min(n_diffraction, len(diff))) + rng.sample(nmr, min(n_nmr, len(nmr))))


def chains_of_entries(provider, entries) -> dict[str, list[str]]:
    wanted = set(entries)
    out: dict[str, list[str]] = defaultdict(list)
    for key in provider.chain_keys():
        if key[:4] in wanted:
            out[key[:4]].append(key[4:])
    return dict(out)


# ---------------------------------------------------------------- full-run candidates
def candidate_chains(provider, resolu: dict[str, float | None], entry_type: dict[str, tuple[str, str]],
                     max_resolution: float) -> pd.DataFrame:
    """Multi-domain chains of protein entries with resolution <= ``max_resolution`` (inclusive)."""
    rows = []
    for key in provider.chain_keys():
        pdb, chain = key[:4], key[4:]
        mol = entry_type.get(pdb, ("", ""))[0]
        res = resolu.get(pdb)
        if mol not in ("prot", "prot-nuc") or res is None or res > max_resolution:
            continue
        rows.append({"pdb_id": pdb, "chain": chain, "resolution": res, "method_class": entry_type[pdb][1]})
    return pd.DataFrame(rows, columns=["pdb_id", "chain", "resolution", "method_class"])


def attach_clusters(cands: pd.DataFrame, entity_info: dict[str, dict], clusters: dict[str, int],
                    n_domains: dict[str, int]) -> pd.DataFrame:
    """Add entity id, UniProt count and 40% cluster id to each candidate chain."""
    ent_col, nunp_col, clu_col, ptype_col = [], [], [], []
    for pdb, chain in zip(cands["pdb_id"], cands["chain"]):
        info = entity_info.get(pdb, {})
        ent = None
        for e in info.get("entities", []):
            if chain in e["auth_asym_ids"]:
                ent = e
                break
        if ent is None and chain == "0":
            prot = [e for e in info.get("entities", []) if e.get("type") == "polypeptide(L)"]
            ent = prot[0] if len(prot) == 1 else None
        ent_col.append(ent["entity"] if ent else None)
        nunp_col.append(len(set(ent["uniprot_ids"])) if ent else 0)
        ptype_col.append(ent.get("type") if ent else None)
        if ent is None:
            clu_col.append(f"nomap:{pdb}{chain}")
        else:
            cid = clusters.get(ent["entity"])
            clu_col.append(f"c{cid}" if cid is not None else f"single:{ent['entity']}")
    out = cands.copy()
    out["entity"] = ent_col
    out["n_uniprot"] = nunp_col
    out["polymer_type"] = ptype_col
    out["cluster"] = clu_col
    out["n_domains"] = [n_domains.get(f"{p}{c}", 0) for p, c in zip(out["pdb_id"], out["chain"])]
    return out


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Order candidates inside each cluster (best first) and add ``rank`` (0-based).

    Exactly one UniProt accession first (likely natural), then best resolution,
    then most domains, then PDB id / chain for determinism. One chain per entity.
    """
    df = df[df["polymer_type"].isna() | (df["polymer_type"] == "polypeptide(L)")].copy()
    df["_unp_bad"] = (df["n_uniprot"] != 1).astype(int)
    df = df.sort_values(["cluster", "_unp_bad", "resolution", "n_domains", "pdb_id", "chain"],
                        ascending=[True, True, True, False, True, True])
    # one chain per entity (homo-oligomer copies are redundant by definition)
    df = df.drop_duplicates(subset=["cluster", "entity", "pdb_id"], keep="first")
    df["rank"] = df.groupby("cluster").cumcount()
    return df.drop(columns="_unp_bad").reset_index(drop=True)
