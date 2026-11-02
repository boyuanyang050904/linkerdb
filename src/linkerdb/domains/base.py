"""Source-independent representation of domain boundaries."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AuthRes:
    """A residue in author (PDB-format) numbering."""

    chain: str
    num: int
    ins: str = ""

    def __str__(self) -> str:
        return f"{self.num}{self.ins}"


@dataclass
class Segment:
    """A contiguous stretch of one domain.

    ``start``/``end`` are SEQRES indices (mmCIF ``label_seq_id``), inclusive.
    Sources that use author numbering fill ``start_auth``/``end_auth`` and leave
    the SEQRES indices to be resolved against the structure.
    """

    start: int | None = None
    end: int | None = None
    start_auth: AuthRes | None = None
    end_auth: AuthRes | None = None

    def resolved(self) -> bool:
        return self.start is not None and self.end is not None


@dataclass
class Domain:
    domain_id: str
    segments: list[Segment]
    classification: str | None = None  # CATH code, SCOPe sccs, Pfam accession
    length: int | None = None          # length reported by the source (CATH: observed residues)


@dataclass
class ChainDomains:
    pdb_id: str
    chain: str                       # author chain id
    source: str                      # e.g. "CATH v4.4.0"
    domains: list[Domain]
    fragments: list[Segment] = field(default_factory=list)
    # alternative SEQRES chopping from a second file of the same source (CATH seqres)
    seqres_domains: list[Domain] | None = None

    @property
    def n_domains(self) -> int:
        return len(self.domains)


def chain_key(pdb_id: str, chain: str) -> str:
    return f"{pdb_id.lower()}{chain}"
