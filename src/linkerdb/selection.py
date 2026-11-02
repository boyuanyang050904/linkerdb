"""Dataset selection: the 50-structure sample, full-run candidates, cluster
representatives and linker-level redundancy removal (docs/design.md §8)."""

from __future__ import annotations

import logging
import random
from collections import defaultdict

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

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


# ---------------------------------------------------------------- linker-level redundancy
def _aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10.0
    al.extend_gap_score = -0.5
    # free end gaps: a 2-domain construct must align inside a 3-domain one without penalty
    names = dir(Align.PairwiseAligner)  # (hasattr would call a getter that can raise)
    for new, old in (("end_insertion_score", "target_end_gap_score"),
                     ("end_deletion_score", "query_end_gap_score")):
        setattr(al, new if new in names else old, 0.0)
    return al


def _clean(seq: str) -> str:
    # letters outside the BLOSUM62 alphabet (U, O, J) are mapped to X
    return "".join(c if c in "ARNDCQEGHILKMFPSTWYVBZX" else "X" for c in seq.upper())


def identity(a: str, b: str, aligner: Align.PairwiseAligner | None = None) -> float:
    """Identical aligned positions / length of the shorter sequence."""
    return align_pair(a, b, aligner)[0]


def align_pair(a: str, b: str, aligner: Align.PairwiseAligner | None = None):
    """Return ``(identity, alignment)`` with identity over the shorter sequence."""
    if not a or not b:
        return 0.0, None
    aligner = aligner or _aligner()
    aln = aligner.align(_clean(a), _clean(b))[0]
    ident = aln.counts().identities
    return ident / min(len(a), len(b)), aln


def map_position(aln, pos: int) -> int | None:
    """Map a 1-based position of the first sequence onto the second (nearest aligned block)."""
    p = pos - 1
    best = None
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        if a0 <= p < a1:
            return b0 + (p - a0) + 1
        dist = a0 - p if p < a0 else p - a1 + 1
        cand = (b0 + 1) if p < a0 else b1
        if best is None or dist < best[0]:
            best = (dist, cand)
    return None if best is None else best[1]


def _overlaps(lo1, hi1, lo2, hi2) -> bool:
    return lo1 <= hi2 and lo2 <= hi1


def _anchors(row) -> tuple[int, int]:
    return int(row["ss_start"]) - 1, int(row["ss_end"]) + 1


def dedup_linkers(linkers: pd.DataFrame, priority: list[str], threshold: float = 0.40,
                  chain_seqs: dict[tuple[str, str], str] | None = None, tolerance: int = 3) -> pd.DataFrame:
    """Remove linkers that are redundant with a better-ranked linker from another chain.

    1. same UniProt accession and overlapping UniProt interval (anchors included);
    2. same (N, C) CATH superfamily pair, full-chain identity >= ``threshold`` over
       the shorter chain (no coverage requirement), and the alignment maps one
       junction onto the other (within ``tolerance`` residues).

    Linkers of the same chain are never compared: they are different junctions.
    Dropped linkers are reported in ``result.attrs['dropped']`` (``{id: reason}``).
    """
    df = linkers.sort_values(priority).reset_index(drop=True)
    dropped: dict[str, str] = {}
    keep_mask = [True] * len(df)
    chain_of = list(zip(df["pdb_id"], df["auth_asym_id"]))
    kept_by_unp: dict[str, list[tuple[int, int, str, tuple]]] = defaultdict(list)
    for i, row in df.iterrows():
        acc, lo, hi = row["uniprot"], row["unp_lo"], row["unp_hi"]
        if acc is None or pd.isna(lo) or pd.isna(hi):
            continue
        lo, hi = int(min(lo, hi)), int(max(lo, hi))
        hit = next((lid for klo, khi, lid, ch in kept_by_unp[acc]
                    if ch != chain_of[i] and _overlaps(lo, hi, klo, khi)), None)
        if hit:
            keep_mask[i] = False
            dropped[row["linker_id"]] = f"same_uniprot_region:{hit}"
        else:
            kept_by_unp[acc].append((lo, hi, row["linker_id"], chain_of[i]))
    if chain_seqs is not None:
        aligner = _aligner()
        cache: dict[tuple, tuple] = {}
        kept_by_pair: dict[tuple, list[int]] = defaultdict(list)
        n_aln = 0
        for i, row in df.iterrows():
            if not keep_mask[i]:
                continue
            key = (row["cath_n"] if isinstance(row["cath_n"], str) else "?",
                   row["cath_c"] if isinstance(row["cath_c"], str) else "?")
            seq_i = chain_seqs.get(chain_of[i])
            hit = None
            if seq_i:
                lo_i, hi_i = _anchors(row)
                for j in kept_by_pair[key]:
                    if chain_of[j] == chain_of[i]:
                        continue
                    seq_j = chain_seqs.get(chain_of[j])
                    if not seq_j:
                        continue
                    pair = (chain_of[i], chain_of[j])
                    if pair not in cache:
                        n_aln += 1
                        cache[pair] = align_pair(seq_i, seq_j, aligner)
                    ident, aln = cache[pair]
                    if ident < threshold or aln is None:
                        continue
                    m_lo, m_hi = map_position(aln, max(lo_i, 1)), map_position(aln, min(hi_i, len(seq_i)))
                    lo_j, hi_j = _anchors(df.loc[j])
                    if None not in (m_lo, m_hi) and _overlaps(m_lo, m_hi, lo_j - tolerance, hi_j + tolerance):
                        hit = f"{df.at[j, 'linker_id']}({ident:.2f})"
                        break
            if hit:
                keep_mask[i] = False
                dropped[row["linker_id"]] = f"homologous_chain:{hit}"
            else:
                kept_by_pair[key].append(i)
        log.info("linker dedup: %d chain alignments", n_aln)
    log.info("linker dedup: %d dropped of %d", len(dropped), len(df))
    out = df[keep_mask].reset_index(drop=True)
    out.attrs["dropped"] = dropped
    return out
