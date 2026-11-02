"""Choose the domain-boundary source: CATH -> SCOPe -> Pfam.

The choice is made once per run (``build_provider``): a source is used only if
its files can be obtained and parsed. ``ChainedProvider`` additionally allows a
per-chain fallback (a chain absent from CATH is looked up in SCOPe, ...); it is
off by default because mixing sources mixes domain definitions.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..config import SCOPE_CLA_URL
from ..download import DownloadError
from .base import ChainDomains
from .cath import CathFormatError, CathProvider, fetch_cath
from .pfam import PfamProvider
from .scope import ScopeProvider

log = logging.getLogger(__name__)

SOURCES = ("cath", "scope", "pfam")


class ProviderUnavailable(RuntimeError):
    pass


def _make(name: str, raw_dir: Path, downloader):
    if name == "cath":
        paths = fetch_cath(raw_dir, downloader)
        prov = CathProvider(paths["cdf"], paths["seqres"], paths["list"],
                            unclassified_path=paths["unclassified"])
        if not prov.chain_keys():
            raise ProviderUnavailable("CATH boundary file contains no multi-domain chains")
        return prov
    if name == "scope":
        path = Path(raw_dir) / "scope" / Path(SCOPE_CLA_URL).name
        downloader.fetch(SCOPE_CLA_URL, path)
        prov = ScopeProvider(path)
        if not prov.chain_keys():
            raise ProviderUnavailable("SCOPe file contains no multi-domain chains")
        return prov
    if name == "pfam":
        return PfamProvider(Path(raw_dir) / "pfam", downloader)
    raise ValueError(f"unknown domain source {name!r}")


def build_provider(raw_dir: Path, downloader, order=SOURCES):
    """Return ``(provider, notes)`` for the first usable source in ``order``."""
    notes: list[str] = []
    for name in order:
        try:
            prov = _make(name, raw_dir, downloader)
        except (DownloadError, CathFormatError, ProviderUnavailable, OSError, ValueError) as exc:
            msg = f"{name} unavailable: {exc}"
            log.warning(msg)
            notes.append(msg)
            continue
        notes.append(f"using {prov.source}")
        return prov, notes
    raise ProviderUnavailable("; ".join(notes))


class ChainedProvider:
    """Per-chain fallback across several providers (first hit wins)."""

    name = "chained"

    def __init__(self, providers) -> None:
        self.providers = list(providers)

    @property
    def source(self) -> str:
        return " > ".join(p.source for p in self.providers)

    def get(self, pdb_id: str, chain: str) -> ChainDomains | None:
        for p in self.providers:
            res = p.get(pdb_id, chain)
            if res is not None:
                return res
        return None
