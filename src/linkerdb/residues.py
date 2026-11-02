"""Residue-name handling: canonical one-letter codes for standard and modified residues."""

from __future__ import annotations

import gemmi

STANDARD_AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# genetically encoded but rare; not counted as "non-standard"
ENCODED_EXTRA = {"SEC": "U", "PYL": "O"}
UNKNOWN = {"UNK", "UNX"}


def one_letter(mon_id: str, parents: dict[str, str] | None = None, _depth: int = 0) -> str:
    """Canonical one-letter code of ``mon_id``.

    Order: standard table -> gemmi's tabulated modified residues (parent amino acid)
    -> ``parents`` (from ``_pdbx_struct_mod_residue``) -> ``'X'``.
    """
    mon_id = mon_id.upper()
    if mon_id in STANDARD_AA:
        return STANDARD_AA[mon_id]
    if mon_id in ENCODED_EXTRA:
        return ENCODED_EXTRA[mon_id]
    if mon_id in UNKNOWN:
        return "X"
    info = gemmi.find_tabulated_residue(mon_id)
    if info is not None and info.is_amino_acid():
        code = info.one_letter_code.upper()
        if code.isalpha() and code != "X":
            return code
    if parents and mon_id in parents and _depth < 3:
        return one_letter(parents[mon_id], parents, _depth + 1)
    return "X"


def is_nonstandard(mon_id: str) -> bool:
    """True for anything that is not one of the 20 standard amino acids (or SEC/PYL)."""
    m = mon_id.upper()
    return m not in STANDARD_AA and m not in ENCODED_EXTRA
