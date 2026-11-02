"""Per-chain structural model built from an mmCIF (or, as fallback, PDB) file.

Everything is indexed by SEQRES position (``label_seq_id``, 1-based). Unmodelled
residues are present with ``observed = False`` and NaN coordinates. Only the
first model is used for coordinates (NMR ensembles), and only the first
alternative conformation is kept.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

from .residues import is_nonstandard, one_letter

BACKBONE = ("N", "CA", "C", "O")


class StructureError(ValueError):
    pass


@dataclass
class Residue:
    seq_id: int
    mon_id: str
    one: str
    auth_num: int | None
    ins: str
    observed: bool

    @property
    def nonstandard(self) -> bool:
        return is_nonstandard(self.mon_id)

    @property
    def auth(self) -> str:
        return "" if self.auth_num is None else f"{self.auth_num}{self.ins}"


@dataclass
class RefSegment:
    """One ``_struct_ref_seq`` row: SEQRES ``seq_beg..seq_end`` aligned to ``db_beg..db_end``."""

    db_name: str
    accession: str
    seq_beg: int
    seq_end: int
    db_beg: int
    db_end: int


@dataclass
class SeqDif:
    seq_id: int
    mon_id: str
    db_mon_id: str
    details: str


@dataclass
class EntryInfo:
    pdb_id: str
    method: str
    resolution: float | None
    n_models: int
    fmt: str = "cif"


@dataclass
class Chain:
    entry: EntryInfo
    auth_chain: str
    label_asym: str
    entity_id: str
    polymer_type: str
    description: str
    taxid: str
    residues: list[Residue]
    coords: np.ndarray  # (N, 4, 3) N/CA/C/O, NaN where missing
    refs: list[RefSegment] = field(default_factory=list)
    seq_difs: list[SeqDif] = field(default_factory=list)
    annotated_ss: str = ""  # from _struct_conf/_struct_sheet_range ('H', 'E', '-')
    microheterogeneity: bool = False
    ref_source: str = "mmcif"  # 'sifts' when UniProt segments were replaced from SIFTS

    def use_sifts(self, segments: list[RefSegment] | None) -> None:
        """Replace UniProt reference segments by SIFTS ones (non-UniProt refs are kept)."""
        if not segments:
            return
        self.refs = [r for r in self.refs if r.db_name != "UNP"] + list(segments)
        self.ref_source = "sifts"

    # derived ------------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.residues)

    @property
    def sequence(self) -> str:
        return "".join(r.one for r in self.residues)

    @property
    def observed(self) -> np.ndarray:
        return np.array([r.observed for r in self.residues], dtype=bool)

    def res(self, seq_id: int) -> Residue:
        return self.residues[seq_id - 1]

    def auth_index(self) -> dict[tuple[int, str], int]:
        """``(auth_num, ins) -> seq_id`` (first occurrence wins)."""
        out: dict[tuple[int, str], int] = {}
        for r in self.residues:
            if r.auth_num is not None:
                out.setdefault((r.auth_num, r.ins), r.seq_id)
        return out

    def uniprot_accessions(self) -> list[str]:
        return sorted({r.accession for r in self.refs if r.db_name == "UNP"})

    def unp_position(self, seq_id: int) -> tuple[str, int] | None:
        for ref in self.refs:
            if ref.db_name == "UNP" and ref.seq_beg <= seq_id <= ref.seq_end:
                return ref.accession, ref.db_beg + (seq_id - ref.seq_beg)
        return None

    def breaks(self, cn_max: float = 2.0, ca_max: float = 4.2) -> np.ndarray:
        """Boolean array ``b`` of length N-1: ``b[i]`` is True if the backbone between
        SEQRES positions ``i+1`` and ``i+2`` is broken although both are observed."""
        c = self.coords
        if self.n < 2:
            return np.zeros(0, dtype=bool)
        cn = np.linalg.norm(c[:-1, 2] - c[1:, 0], axis=1)
        caca = np.linalg.norm(c[:-1, 1] - c[1:, 1], axis=1)
        obs = self.observed
        both = obs[:-1] & obs[1:]
        with np.errstate(invalid="ignore"):
            brk = np.where(np.isnan(cn), caca > ca_max, cn > cn_max)
        # both observed but CA missing too (e.g. only side chain atoms): treat as break
        brk = np.where(np.isnan(cn) & np.isnan(caca), True, brk)
        return both & brk


# ---------------------------------------------------------------- parsing helpers
def _table(block: gemmi.cif.Block, prefix: str, tags: list[str]) -> list[list[str | None]]:
    tab = block.find(prefix, tags)
    rows = []
    for row in tab:
        vals = []
        for i in range(len(tags)):
            if not row.has(i) or gemmi.cif.is_null(row[i]):
                vals.append(None)
            else:
                vals.append(row.str(i))
        rows.append(vals)
    return rows


def _int(v: str | None) -> int | None:
    try:
        return int(v) if v is not None else None
    except ValueError:
        return None


def _entry_info(block: gemmi.cif.Block, st: gemmi.Structure, pdb_id: str, fmt: str) -> EntryInfo:
    methods = [r[0] for r in _table(block, "_exptl.", ["method"]) if r[0]] if block is not None else []
    if not methods and "_exptl.method" in st.info:
        methods = [st.info["_exptl.method"]]
    method = ";".join(methods)
    res = st.resolution if st.resolution and st.resolution > 0 else None
    return EntryInfo(pdb_id.lower(), method.upper(), res, len(st), fmt)


def _backbone(st: gemmi.Structure, auth_chain: str, n: int) -> tuple[np.ndarray, dict[int, str]]:
    coords = np.full((n, 4, 3), np.nan)
    names: dict[int, str] = {}
    model = st[0]
    for ch in model:
        if ch.name != auth_chain:
            continue
        for res in ch:
            if res.entity_type != gemmi.EntityType.Polymer or res.label_seq is None:
                continue
            i = res.label_seq - 1
            if not 0 <= i < n or i in names:
                continue
            names[i] = res.name
            for k, an in enumerate(BACKBONE):
                atom = res.find_atom(an, "*")
                if atom is not None:
                    p = atom.pos
                    coords[i, k] = (p.x, p.y, p.z)
    return coords, names


def _annotated_ss(block: gemmi.cif.Block, label_asym: str, n: int) -> str:
    ss = ["-"] * n
    for typ, asym, b, e in _table(block, "_struct_conf.", ["conf_type_id", "beg_label_asym_id",
                                                           "beg_label_seq_id", "end_label_seq_id"]):
        if asym == label_asym and typ and typ.startswith("HELX") and _int(b) and _int(e):
            for i in range(_int(b) - 1, min(_int(e), n)):
                ss[i] = "H"
    for asym, b, e in _table(block, "_struct_sheet_range.", ["beg_label_asym_id", "beg_label_seq_id",
                                                             "end_label_seq_id"]):
        if asym == label_asym and _int(b) and _int(e):
            for i in range(_int(b) - 1, min(_int(e), n)):
                ss[i] = "E"
    return "".join(ss)


# ---------------------------------------------------------------- mmCIF
def load_mmcif(path: Path | str, pdb_id: str | None = None) -> "EntryModel":
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    pdb_id = (pdb_id or block.name).lower()
    st = gemmi.make_structure_from_block(block)
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    info = _entry_info(block, st, pdb_id, "cif")

    entity_type = {r[0]: r[1] for r in _table(block, "_entity_poly.", ["entity_id", "type"])}
    entity_can = {r[0]: (r[1] or "").replace("\n", "").replace(" ", "")
                  for r in _table(block, "_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code_can"])}
    entity_desc = {r[0]: r[1] or "" for r in _table(block, "_entity.", ["id", "pdbx_description"])}
    taxid: dict[str, str] = {}
    for prefix, tag in (("_entity_src_gen.", "pdbx_gene_src_ncbi_taxonomy_id"),
                        ("_entity_src_nat.", "pdbx_ncbi_taxonomy_id"),
                        ("_pdbx_entity_src_syn.", "ncbi_taxonomy_id")):
        for eid, tx in _table(block, prefix, ["entity_id", tag]):
            if eid and tx:
                taxid.setdefault(eid, tx)
    parents = {r[0]: r[1] for r in _table(block, "_pdbx_struct_mod_residue.",
                                          ["label_comp_id", "parent_comp_id"]) if r[0] and r[1]}

    refs_db = {r[0]: (r[1] or "", r[2] or "") for r in _table(block, "_struct_ref.",
                                                              ["id", "db_name", "pdbx_db_accession"])}
    ref_rows = _table(block, "_struct_ref_seq.", ["ref_id", "pdbx_strand_id", "seq_align_beg", "seq_align_end",
                                                  "pdbx_db_accession", "db_align_beg", "db_align_end"])
    dif_rows = _table(block, "_struct_ref_seq_dif.", ["pdbx_pdb_strand_id", "seq_num", "mon_id", "db_mon_id",
                                                      "details"])

    scheme = _table(block, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "pdb_seq_num",
                                                      "pdb_ins_code", "pdb_strand_id", "auth_mon_id", "hetero"])
    if not scheme:
        raise StructureError(f"{pdb_id}: no _pdbx_poly_seq_scheme (not a polymer entry?)")
    per_chain: dict[str, list[list[str | None]]] = {}
    for row in scheme:
        per_chain.setdefault(row[6], []).append(row)

    chains: dict[str, Chain] = {}
    for auth_chain, rows in per_chain.items():
        label_asym, entity_id = rows[0][0], rows[0][1]
        n = max(_int(r[2]) or 0 for r in rows)
        seen: dict[int, list[str | None]] = {}
        micro = False
        for r in rows:
            sid = _int(r[2])
            if sid in seen:
                micro = True
                continue
            seen[sid] = r
        coords, model_names = _backbone(st, auth_chain, n)
        can = entity_can.get(entity_id, "")
        use_can = len(can) == n and can.isalpha()
        residues = []
        for sid in range(1, n + 1):
            r = seen.get(sid)
            if r is None:  # should not happen; keep the position as unknown
                residues.append(Residue(sid, "UNK", "X", None, "", False))
                continue
            mon = r[3] or "UNK"
            # if the modelled residue is a different microheterogeneity variant, prefer it
            mon = model_names.get(sid - 1, mon)
            one = can[sid - 1].upper() if use_can and mon == r[3] else one_letter(mon, parents)
            # "observed" = present in the first model (the coordinate source)
            observed = (sid - 1) in model_names
            residues.append(Residue(sid, mon, one, _int(r[4]), r[5] or "", observed))
        refs = []
        for ref_id, strand, sb, se, acc, db, de in ref_rows:
            if strand != auth_chain or None in (sb, se, db, de):
                continue
            db_name, ref_acc = refs_db.get(ref_id, ("", ""))
            refs.append(RefSegment(db_name, acc or ref_acc, int(sb), int(se), int(db), int(de)))
        difs = [SeqDif(int(sn), mon or "", dbm or "", (det or "").lower())
                for strand, sn, mon, dbm, det in dif_rows if strand == auth_chain and _int(sn) is not None]
        chains[auth_chain] = Chain(
            entry=info, auth_chain=auth_chain, label_asym=label_asym or "", entity_id=entity_id or "",
            polymer_type=entity_type.get(entity_id, ""), description=entity_desc.get(entity_id, ""),
            taxid=taxid.get(entity_id, ""), residues=residues, coords=coords, refs=refs, seq_difs=difs,
            annotated_ss=_annotated_ss(block, label_asym, n), microheterogeneity=micro)
    return EntryModel(info, chains)


# ---------------------------------------------------------------- PDB format fallback
def load_pdb(path: Path | str, pdb_id: str) -> "EntryModel":
    """PDB-format fallback: SEQRES is aligned to the model by gemmi.

    Author numbers of unmodelled residues are unknown (``None``); there is no
    ``_struct_ref_seq`` information beyond DBREF, which gemmi does not expose,
    so natural-protein checks will fail closed (``no_uniprot``).
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.assign_label_seq_id(force=True)
    info = _entry_info(None, st, pdb_id, "pdb")
    chains: dict[str, Chain] = {}
    for ch in st[0]:
        poly = ch.get_polymer()
        if len(poly) == 0 or ch.name in chains:
            continue
        ent = st.get_entity_of(poly)
        seq = list(ent.full_sequence) if ent is not None and ent.full_sequence else [r.name for r in poly]
        seq = [gemmi.Entity.first_mon(s) for s in seq]
        n = len(seq)
        coords, names = _backbone(st, ch.name, n)
        auth = {r.label_seq: r for r in poly if r.label_seq is not None}
        residues = []
        for sid in range(1, n + 1):
            r = auth.get(sid)
            mon = r.name if r is not None else seq[sid - 1]
            residues.append(Residue(sid, mon, one_letter(mon), r.seqid.num if r is not None else None,
                                    (r.seqid.icode.strip() if r is not None else ""), r is not None))
        ptype = "polypeptide(L)" if ent is None or ent.polymer_type == gemmi.PolymerType.PeptideL else str(
            ent.polymer_type)
        chains[ch.name] = Chain(entry=info, auth_chain=ch.name, label_asym=poly.subchain_id(),
                                entity_id=ent.name if ent is not None else "", polymer_type=ptype,
                                description="", taxid="", residues=residues, coords=coords)
    return EntryModel(info, chains)


@dataclass
class EntryModel:
    info: EntryInfo
    chains: dict[str, Chain]


def load_structure(path: Path | str, pdb_id: str, fmt: str = "cif") -> EntryModel:
    return load_mmcif(path, pdb_id) if fmt == "cif" else load_pdb(path, pdb_id)
