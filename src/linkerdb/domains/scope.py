"""SCOPe domain boundaries (fallback 1).

``dir.cla.scope.<ver>.txt`` lines look like::

    d1cnsa1  1cns  A:1-87,A:146-243  d.2.1.1  ...
    d1ux8a_  1ux8  A:                 a.1.1.1  ...   (whole chain)

Ranges use author numbering, optionally with insertion codes (``A:12A-140``)
and negative numbers (``A:-3-120``).
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from ..config import SCOPE_VERSION
from .base import AuthRes, ChainDomains, Domain, Segment

_RANGE = re.compile(r"^(?P<s>-?\d+)(?P<si>[A-Za-z]?)-(?P<e>-?\d+)(?P<ei>[A-Za-z]?)$")


def parse_scope_region(region: str) -> list[tuple[str, Segment | None]]:
    """``"A:1-87,A:146-243"`` -> ``[("A", Segment), ...]``; ``"A:"`` -> ``[("A", None)]``."""
    out: list[tuple[str, Segment | None]] = []
    for part in region.split(","):
        if ":" not in part:
            out.append(("", None))
            continue
        chain, rng = part.split(":", 1)
        if not rng:
            out.append((chain, None))
            continue
        m = _RANGE.match(rng)
        if not m:
            raise ValueError(f"unparsable SCOPe range {part!r}")
        out.append((chain, Segment(start_auth=AuthRes(chain, int(m["s"]), m["si"]),
                                   end_auth=AuthRes(chain, int(m["e"]), m["ei"]))))
    return out


class ScopeProvider:
    name = "SCOPe"

    def __init__(self, cla_path: Path, version: str = SCOPE_VERSION, min_domains: int = 2) -> None:
        self.version = version
        by_chain: dict[str, list[Domain]] = defaultdict(list)
        with open(cla_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                t = line.split("\t") if "\t" in line else line.split()
                if len(t) < 4:
                    continue
                sid, pdb, region, sccs = t[0], t[1], t[2], t[3]
                try:
                    parts = parse_scope_region(region)
                except ValueError:
                    continue
                chains = {c for c, _ in parts}
                if len(chains) != 1 or any(seg is None for _, seg in parts):
                    continue  # whole-chain or multi-chain domain: no internal boundary
                chain = chains.pop()
                by_chain[f"{pdb.lower()}{chain}"].append(Domain(sid, [s for _, s in parts], sccs))
        self._chains = {k: v for k, v in by_chain.items() if len(v) >= min_domains}

    @property
    def source(self) -> str:
        return f"SCOPe {self.version}"

    def chain_keys(self) -> list[str]:
        return sorted(self._chains)

    def entries(self) -> set[str]:
        return {k[:4] for k in self._chains}

    def get(self, pdb_id: str, chain: str) -> ChainDomains | None:
        doms = self._chains.get(f"{pdb_id.lower()}{chain}")
        if doms is None:
            return None
        return ChainDomains(pdb_id.lower(), chain, self.source, sorted(doms, key=lambda d: d.domain_id))
