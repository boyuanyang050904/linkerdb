"""CATH domain boundaries.

Two CATH files describe the same chopping in different numbering schemes:

* ``cath-domain-boundaries.txt`` - CATH Domall Format 2.0 (CDF) in author
  numbering with insertion codes, including unassigned *fragments*;
* ``cath-domain-boundaries-seqreschopping.txt`` - SEQRES indices.

Both are parsed so that the pipeline can cross-validate the mapping onto the
current mmCIF file (docs/design.md §7.4).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

from ..config import CATH_BASE, CATH_FILES, CATH_VERSION
from .base import AuthRes, ChainDomains, Domain, Segment

log = logging.getLogger(__name__)


class CathFormatError(ValueError):
    pass


def _ins(tok: str) -> str:
    return "" if tok == "-" else tok


def parse_cdf_line(line: str) -> tuple[str, list[list[Segment]], list[Segment]]:
    """Parse one CDF 2.0 line into ``(chain_name, domains, fragments)``.

    ``domains`` is a list (one entry per domain, in file order) of segment lists.
    """
    tok = line.split()
    if len(tok) < 3 or not tok[1].startswith("D") or not tok[2].startswith("F"):
        raise CathFormatError(f"not a CDF line: {line!r}")
    name = tok[0]
    try:
        n_dom = int(tok[1][1:])
        n_frag = int(tok[2][1:])
    except ValueError as exc:
        raise CathFormatError(f"bad domain/fragment count: {line!r}") from exc
    i = 3
    domains: list[list[Segment]] = []
    try:
        for _ in range(n_dom):
            n_seg = int(tok[i])
            i += 1
            segs = []
            for _ in range(n_seg):
                c1, s, si, c2, e, ei = tok[i:i + 6]
                i += 6
                segs.append(Segment(start_auth=AuthRes(c1, int(s), _ins(si)),
                                    end_auth=AuthRes(c2, int(e), _ins(ei))))
            domains.append(segs)
        fragments = []
        for _ in range(n_frag):
            c1, s, si, c2, e, ei, _nr = tok[i:i + 7]
            i += 7
            fragments.append(Segment(start_auth=AuthRes(c1, int(s), _ins(si)),
                                     end_auth=AuthRes(c2, int(e), _ins(ei))))
    except (ValueError, IndexError) as exc:
        raise CathFormatError(f"truncated or malformed CDF line: {line!r}") from exc
    if i != len(tok):
        raise CathFormatError(f"trailing tokens in CDF line: {line!r}")
    return name, domains, fragments


def parse_seqres_ranges(ranges: str) -> list[Segment]:
    segs = []
    for part in ranges.split(","):
        a, b = part.rsplit("-", 1) if part.count("-") == 1 else _split_signed(part)
        segs.append(Segment(start=int(a), end=int(b)))
    return segs


def _split_signed(part: str) -> tuple[str, str]:
    # tolerate negative numbers such as "-3-10" (never seen in SEQRES files, but cheap)
    idx = part.index("-", 1)
    return part[:idx], part[idx + 1:]


class CathProvider:
    name = "CATH"

    def __init__(self, cdf_path: Path, seqres_path: Path | None = None, list_path: Path | None = None,
                 version: str = CATH_VERSION, min_domains: int = 2, unclassified_path: Path | None = None) -> None:
        self.version = version
        self.min_domains = min_domains
        self._cdf: dict[str, str] = {}
        self._seqres: dict[str, list[tuple[str, str]]] = defaultdict(list)
        self._info: dict[str, tuple[str, int]] = {}
        self._load_cdf(Path(cdf_path))
        if seqres_path is not None and Path(seqres_path).exists():
            self._load_seqres(Path(seqres_path))
        if list_path is not None and Path(list_path).exists():
            self._load_list(Path(list_path))
        if unclassified_path is not None and Path(unclassified_path).exists():
            # chopped but not yet classified domains: keep the length, no CATH code
            self._load_list(Path(unclassified_path), classified=False)

    @property
    def source(self) -> str:
        return f"CATH {self.version}"

    # loading -------------------------------------------------------------
    def _load_cdf(self, path: Path) -> None:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                name = line[:5]
                head = line.split(None, 2)
                if len(head) < 2 or not head[1].startswith("D"):
                    continue
                if int(head[1][1:]) >= self.min_domains:
                    self._cdf[name] = line.rstrip("\n")

    def _load_seqres(self, path: Path) -> None:
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) != 2:
                    continue
                did, ranges = parts
                if did[:5] in self._cdf:
                    self._seqres[did[:5]].append((did, ranges))

    def _load_list(self, path: Path, classified: bool = True) -> None:
        # CLF columns are fixed width and large counts can merge with their neighbour
        # (e.g. "1110732"), so only the leading C.A.T.H and the trailing length /
        # resolution columns are read by position.
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                t = line.split()
                if len(t) < 8 or t[0][:5] not in self._cdf:
                    continue
                try:
                    length = int(t[-2])
                except ValueError:
                    continue
                code = ".".join(t[1:5]) if classified else None
                self._info.setdefault(t[0], (code, length))

    # queries -------------------------------------------------------------
    def chain_keys(self) -> list[str]:
        """All chains with >= ``min_domains`` domains, as ``pdbid+chain`` (e.g. ``1cnsA``)."""
        return sorted(self._cdf)

    def entries(self) -> set[str]:
        return {k[:4] for k in self._cdf}

    def get(self, pdb_id: str, chain: str) -> ChainDomains | None:
        key = f"{pdb_id.lower()}{chain}"
        line = self._cdf.get(key)
        if line is None:
            return None
        _, dom_segs, frags = parse_cdf_line(line)
        domains = []
        for k, segs in enumerate(dom_segs, start=1):
            did = f"{key}{k:02d}"
            code, length = self._info.get(did, (None, None))
            domains.append(Domain(did, segs, code, length))
        seqres_domains = None
        if key in self._seqres:
            seqres_domains = []
            for did, ranges in sorted(self._seqres[key]):
                code, length = self._info.get(did, (None, None))
                seqres_domains.append(Domain(did, parse_seqres_ranges(ranges), code, length))
        return ChainDomains(pdb_id.lower(), chain, self.source, domains, frags, seqres_domains)


def cath_paths(raw_dir: Path) -> dict[str, Path]:
    return {k: Path(raw_dir) / "cath" / v for k, v in CATH_FILES.items()}


def fetch_cath(raw_dir: Path, downloader) -> dict[str, Path]:
    paths = cath_paths(raw_dir)
    for key, path in paths.items():
        downloader.fetch(f"{CATH_BASE}/{CATH_FILES[key]}", path)
    return paths
