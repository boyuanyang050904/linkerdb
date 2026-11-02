"""SIFTS chain-level UniProt mapping (``pdb_chain_uniprot.tsv.gz``).

SIFTS is curated by PDBe and updated weekly, so it also covers entries whose
mmCIF ``_struct_ref`` points to a non-UniProt database (e.g. ``PDB`` self
references). ``RES_BEG``/``RES_END`` are SEQRES indices; ``SP_BEG``/``SP_END``
UniProt positions. It is the primary UniProt mapping; ``_struct_ref_seq`` from
the mmCIF file is the fallback for chains SIFTS does not list.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

from .config import SIFTS_UNIPROT_URL
from .structure import RefSegment


def sifts_path(raw_dir: Path) -> Path:
    return Path(raw_dir) / "sifts" / "pdb_chain_uniprot.tsv.gz"


def fetch_sifts(raw_dir: Path, downloader) -> Path:
    path = sifts_path(raw_dir)
    downloader.fetch(SIFTS_UNIPROT_URL, path)
    return path


def load_sifts(path: Path, entries: set[str] | None = None) -> dict[str, dict[str, list[RefSegment]]]:
    """``{pdb_id: {auth_chain: [RefSegment, ...]}}``, optionally restricted to ``entries``."""
    out: dict[str, dict[str, list[RefSegment]]] = defaultdict(lambda: defaultdict(list))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("PDB\t"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            pdb = t[0].lower()
            if entries is not None and pdb not in entries:
                continue
            try:
                seg = RefSegment("UNP", t[2], int(t[3]), int(t[4]), int(t[7]), int(t[8]))
            except ValueError:
                continue
            out[pdb][t[1]].append(seg)
    return {p: dict(c) for p, c in out.items()}
