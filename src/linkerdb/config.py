"""Central configuration: data locations, remote URLs and rule parameters.

Everything that influences the content of the output dataset lives here so a run
can be reproduced from the values recorded in ``run_info.json``.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path

CATH_VERSION = "v4.4.0"
CATH_BASE = "https://download.cathdb.info/cath/releases/latest-release/cath-classification-data"
CATH_FILES = {
    "cdf": "cath-domain-boundaries.txt",
    "seqres": "cath-domain-boundaries-seqreschopping.txt",
    "list": "cath-domain-list.txt",
    "unclassified": "cath-unclassified-list.txt",
}

SCOPE_VERSION = "2.08"
SCOPE_CLA_URL = "https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt"

PDBE_PFAM_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb}"
SIFTS_UNIPROT_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"

WWPDB_INDEX = "https://files.wwpdb.org/pub/pdb/derived_data"
RESOLU_URL = f"{WWPDB_INDEX}/index/resolu.idx"
ENTRY_TYPE_URL = f"{WWPDB_INDEX}/pdb_entry_type.txt"
RCSB_CLUSTERS_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{pct}.txt"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

# Structure mirrors, tried in order. ``{id}`` is upper-case, ``{lid}`` lower-case,
# ``{mid}`` the two middle characters used by the divided layout.
MMCIF_MIRRORS = (
    "https://files.rcsb.org/download/{id}.cif.gz",
    "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/{mid}/{lid}.cif.gz",
    "https://www.ebi.ac.uk/pdbe/entry-files/download/{lid}.cif",
)
PDB_FORMAT_MIRRORS = (
    "https://files.rcsb.org/download/{id}.pdb.gz",
)


def default_data_dir() -> Path:
    return Path(os.environ.get("LINKERDB_DATA", Path.cwd() / "data"))


@dataclass(frozen=True)
class Rules:
    """Parameters of the linker definition and QC (see docs/design.md §4–§7)."""

    max_ext: int = 15              # max SS-based extension per side (sample: P95 of uncapped extension = 15)
    flank: int = 10                # flanking residues recorded on each side
    long_gap: int = 30             # CATH gap longer than this -> 'long_gap'
    peptide_cn_max: float = 2.0    # C(i)-N(i+1) above this = chain break (Å)
    ca_ca_max: float = 4.2         # fallback CA-CA break threshold (Å)
    ca_step_max: float = 3.80      # max CA-CA distance per sequence step (geometry QC) ...
    ca_tol: float = 0.50           # ... plus this tolerance (real trans CA-CA spans 3.75-3.95 Å)
    domain_len_tol: int = 2        # tolerated |observed - CATH length| per domain
    context_window: int = 20       # flanking residues stored in window_seq (context for users)
    max_resolution: float = 2.5    # inclusive
    tries_per_cluster: int = 3

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class Settings:
    data_dir: Path = field(default_factory=default_data_dir)
    rules: Rules = field(default_factory=Rules)
    timeout: float = 60.0
    retries: int = 4
    workers: int = 8

    @property
    def raw_dir(self) -> Path:
        return Path(self.data_dir) / "raw"

    @property
    def work_dir(self) -> Path:
        return Path(self.data_dir) / "work"

    def structure_path(self, pdb_id: str, fmt: str = "cif") -> Path:
        lid = pdb_id.lower()
        ext = "cif.gz" if fmt == "cif" else "pdb.gz"
        return self.raw_dir / ("mmcif" if fmt == "cif" else "pdb") / lid[1:3] / f"{lid}.{ext}"
