"""Reference synthetic linkers used in fusion-protein engineering, for comparison
against the natural linkers in the database.

Only two families are in common use:

* **(G4S)n** — ``GGGGS`` repeats: the default *flexible* linker (random coil).
* **(EAAAK)n** — ``EAAAK`` repeats: the common *rigid* linker (a single alpha-helix).

Both were designed to be inert, so neither carries a sequence liability. Their
"reach" (how far apart they can hold two domains) is modelled simply: an EAAAK
helix extends ~1.5 Å per residue along its axis (helix rise); (G4S)n is a
flexible coil with no fixed length — its end-to-end distance is a broad
distribution, so a single reach number is not meaningful.
"""

from __future__ import annotations

import pandas as pd

from . import developability

HELIX_RISE = 1.5   # Å per residue along an alpha-helix axis
EXTENDED_COIL = 3.5  # Å per residue, fully-extended coil (an upper bound on flexible reach)

ARTIFICIAL = [
    ("(G4S)1", "GGGGS", "flexible"),
    ("(G4S)2", "GGGGSGGGGS", "flexible"),
    ("(G4S)3", "GGGGSGGGGSGGGGS", "flexible"),
    ("(G4S)4", "GGGGSGGGGSGGGGSGGGGS", "flexible"),
    ("(EAAAK)1", "EAAAK", "rigid_helical"),
    ("(EAAAK)2", "EAAAKEAAAK", "rigid_helical"),
    ("(EAAAK)3", "EAAAKEAAAKEAAAK", "rigid_helical"),
    ("(EAAAK)4", "EAAAKEAAAKEAAAKEAAAK", "rigid_helical"),
]


def artificial_table() -> pd.DataFrame:
    """One row per reference artificial linker with length, class, modelled reach and liabilities."""
    rows = []
    for name, seq, cls in ARTIFICIAL:
        liab = developability.scan(seq)
        rows.append({
            "name": name,
            "sequence": seq,
            "length": len(seq),
            "class": cls,
            "axial_reach_A": round(HELIX_RISE * len(seq), 1) if cls.startswith("rigid") else None,
            "liability_score": developability.liability_score(liab),
            "liability_tags": developability.summary_flags(liab),
        })
    return pd.DataFrame(rows)
