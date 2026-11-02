"""Shared fixtures. Tests run fully offline: any request to a non-local host fails."""

from __future__ import annotations

import shutil
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pytest
import requests

from linkerdb.config import Rules, Settings
from linkerdb.structure import Chain, EntryInfo, RefSegment, Residue, SeqDif

DATA = Path(__file__).parent / "data"
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
         "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
         "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL", "X": "UNK"}


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    real = requests.Session.request

    def guarded(self, method, url, *a, **kw):
        host = urlparse(url).hostname
        if host not in ("127.0.0.1", "localhost"):
            raise RuntimeError(f"network access attempted in tests: {url}")
        return real(self, method, url, *a, **kw)

    monkeypatch.setattr(requests.Session, "request", guarded)


@pytest.fixture
def data_dir(tmp_path) -> Path:
    """A private copy of the fixture cache (CATH/SCOPe/SIFTS/index files + mmCIF)."""
    dst = tmp_path / "data"
    shutil.copytree(DATA / "raw", dst / "raw")
    return dst


@pytest.fixture
def settings(data_dir) -> Settings:
    return Settings(data_dir=data_dir, workers=1, retries=0)


def backbone(n: int, shift_after: int | None = None, shift: float = 10.0) -> np.ndarray:
    """Extended chain with ideal spacing: CA-CA 3.8 Å, C(i)-N(i+1) 1.35 Å."""
    c = np.zeros((n, 4, 3))
    for i in range(n):
        x = 3.8 * i + (shift if shift_after is not None and i >= shift_after else 0.0)
        c[i, 0] = (x, 0, 0)
        c[i, 1] = (x + 1.2, 0.5, 0)
        c[i, 2] = (x + 2.45, 0, 0)
        c[i, 3] = (x + 2.45, 1.2, 0)
    return c


def make_chain(seq: str, observed: str | None = None, mon: dict[int, str] | None = None,
               shift_after: int | None = None, refs=None, difs=None, auth_offset: int = 0,
               polymer_type: str = "polypeptide(L)", pdb_id: str = "9xyz", chain: str = "A") -> Chain:
    """Synthetic chain; ``observed`` is a string of '1'/'0'; ``mon`` overrides 3-letter codes (1-based)."""
    n = len(seq)
    obs = [c == "1" for c in (observed or "1" * n)]
    mon = mon or {}
    residues = []
    for i, aa in enumerate(seq, start=1):
        residues.append(Residue(i, mon.get(i, THREE[aa]), aa, i + auth_offset, "", obs[i - 1]))
    coords = backbone(n, shift_after)
    coords[~np.array(obs)] = np.nan
    if refs is None:
        refs = [RefSegment("UNP", "P12345", 1, n, 101, 100 + n)]
    return Chain(EntryInfo(pdb_id, "X-RAY DIFFRACTION", 2.0, 1), chain, chain, "1", polymer_type, "test protein",
                 "9606", residues, coords, refs, difs or [])


@pytest.fixture
def rules() -> Rules:
    return Rules()


__all__ = ["make_chain", "backbone", "DATA", "SeqDif", "RefSegment"]
