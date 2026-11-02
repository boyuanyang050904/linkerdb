"""Pfam-based (sequence) domain boundaries via the PDBe/SIFTS API (fallback 2).

Pfam regions are sequence-family envelopes rather than structural domains, so
regions may overlap or nest (e.g. the 7tm_1 envelope of 2RH1 spans the inserted
T4 lysozyme). Overlaps are resolved as follows:

* a region fully contained in another one is treated as an *insertion*: the
  container is split around it (it becomes a discontinuous domain);
* partially overlapping regions are cut at the midpoint of the overlap.

PDBe ``residue_number`` is the SEQRES index (``label_seq_id``).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from ..config import PDBE_PFAM_URL
from ..download import NotFoundError
from .base import ChainDomains, Domain, Segment


def _resolve_overlaps(regions: list[tuple[int, int, str, str]]) -> list[Domain]:
    """``regions``: ``(start, end, accession, name)`` on one chain."""
    regions = sorted(regions, key=lambda r: (r[0], -r[1]))
    # drop exact duplicates
    uniq = []
    for r in regions:
        if not uniq or (r[0], r[1]) != (uniq[-1][0], uniq[-1][1]):
            uniq.append(r)
    segs: list[list[list[int]]] = [[[s, e]] for s, e, _, _ in uniq]
    for i, (si, ei, _, _) in enumerate(uniq):
        for j, (sj, ej, _, _) in enumerate(uniq):
            if i == j:
                continue
            if si <= sj and ej <= ei and (si, ei) != (sj, ej):
                # j nested in i -> cut i around j
                new = []
                for a, b in segs[i]:
                    if b < sj or a > ej:
                        new.append([a, b])
                        continue
                    if a < sj:
                        new.append([a, sj - 1])
                    if b > ej:
                        new.append([ej + 1, b])
                segs[i] = new
    # partial overlaps between consecutive pieces of different domains
    pieces = sorted((a, b, k) for k, ss in enumerate(segs) for a, b in ss)
    fixed: dict[int, list[list[int]]] = defaultdict(list)
    for idx, (a, b, k) in enumerate(pieces):
        if idx + 1 < len(pieces):
            na, nb, nk = pieces[idx + 1]
            if nk != k and na <= b:
                mid = (na + b) // 2
                b = mid
                pieces[idx + 1] = (mid + 1, nb, nk)
        if a <= b:
            fixed[k].append([a, b])
    out = []
    for k, (s, e, acc, name) in enumerate(uniq):
        if fixed.get(k):
            out.append(Domain(f"{acc}:{s}-{e}", [Segment(start=a, end=b) for a, b in fixed[k]], acc))
    return out


def parse_pdbe_pfam(payload: dict, pdb_id: str) -> dict[str, list[Domain]]:
    """Return ``{auth_chain: [Domain, ...]}`` from a PDBe ``mappings/pfam`` response."""
    entry = payload.get(pdb_id.lower(), {}).get("Pfam", {})
    regions: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    for acc, rec in entry.items():
        for m in rec.get("mappings", []):
            regions[m["chain_id"]].append(
                (int(m["start"]["residue_number"]), int(m["end"]["residue_number"]), acc, rec.get("name", "")))
    return {ch: _resolve_overlaps(r) for ch, r in regions.items()}


class PfamProvider:
    name = "Pfam"

    def __init__(self, cache_dir: Path, downloader, min_domains: int = 2) -> None:
        self.cache_dir = Path(cache_dir)
        self.downloader = downloader
        self.min_domains = min_domains

    @property
    def source(self) -> str:
        return "Pfam (PDBe SIFTS mapping)"

    def _payload(self, pdb_id: str) -> dict:
        path = self.cache_dir / f"{pdb_id.lower()}.json"
        if path.exists():
            return json.loads(path.read_text())
        try:
            data = self.downloader.get_json(PDBE_PFAM_URL.format(pdb=pdb_id.lower()))
        except NotFoundError:
            data = {}
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))
        return data

    def chain_keys(self) -> list[str]:  # pragma: no cover - not enumerable offline
        raise NotImplementedError("Pfam provider is queried per entry; it cannot enumerate chains")

    def get(self, pdb_id: str, chain: str) -> ChainDomains | None:
        doms = parse_pdbe_pfam(self._payload(pdb_id), pdb_id).get(chain)
        if not doms or len(doms) < self.min_domains:
            return None
        return ChainDomains(pdb_id.lower(), chain, self.source, doms)
