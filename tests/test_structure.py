import json

import gemmi
import numpy as np
import pytest

from linkerdb import dssp
from linkerdb.residues import is_nonstandard, one_letter
from linkerdb.structure import RefSegment, StructureError, load_mmcif, load_pdb, load_structure

from .conftest import DATA, backbone, make_chain


def cif(pid):
    return DATA / "raw" / "mmcif" / pid[1:3] / f"{pid}.cif.gz"


# ---------------------------------------------------------------- residues
@pytest.mark.parametrize("mon,expected", [("ALA", "A"), ("MSE", "M"), ("SEP", "S"), ("TPO", "T"), ("PTR", "Y"),
                                          ("CSO", "C"), ("SEC", "U"), ("PYL", "O"), ("UNK", "X"), ("HOH", "X")])
def test_one_letter(mon, expected):
    assert one_letter(mon) == expected


def test_one_letter_parent_table_and_flags():
    assert one_letter("ZZZ") == "X"
    assert one_letter("ZZZ", {"ZZZ": "LYS"}) == "K"
    assert one_letter("ZZZ", {"ZZZ": "ZZZ"}) == "X"  # self-loop terminates
    assert is_nonstandard("MSE") and not is_nonstandard("ALA") and not is_nonstandard("SEC")


# ---------------------------------------------------------------- mmCIF parsing
def test_load_mmcif_basic():
    m = load_mmcif(cif("1cns"), "1cns")
    assert m.info.method == "X-RAY DIFFRACTION" and m.info.resolution == pytest.approx(1.91)
    assert m.info.n_models == 1 and set(m.chains) == {"A", "B"}
    ch = m.chains["A"]
    assert ch.n == 243 and ch.observed.all() and ch.polymer_type == "polypeptide(L)"
    assert ch.uniprot_accessions() == ["P23951"] and ch.description == "CHITINASE"
    assert ch.sequence.startswith("SVSSIVSRAQFDRMLL")
    assert ch.unp_position(1) == ("P23951", ch.refs[0].db_beg)
    assert ch.unp_position(10_000) is None
    assert ch.annotated_ss.count("H") > 50
    assert not ch.breaks().any()


def test_missing_residues_keep_author_numbers():
    ch = load_mmcif(cif("2dch"), "2dch").chains["X"]
    missing = [r for r in ch.residues if not r.observed]
    assert missing and all(r.auth_num is not None for r in missing)
    assert np.isnan(ch.coords[missing[0].seq_id - 1]).all()
    idx = ch.auth_index()
    assert idx[(missing[0].auth_num, missing[0].ins)] == missing[0].seq_id


def test_nmr_uses_first_model():
    m = load_mmcif(cif("2j4n"), "2j4n")
    assert m.info.n_models == 6 and "NMR" in m.info.method and m.info.resolution is None
    ch = next(iter(m.chains.values()))
    st = gemmi.read_structure(str(cif("2j4n")))
    first = st[0][ch.auth_chain].get_polymer()
    ca = next(a for a in first[0] if a.name == "CA").pos
    assert np.allclose(ch.coords[first[0].label_seq - 1, 1], (ca.x, ca.y, ca.z))


def test_chimera_refs_and_seq_difs():
    ch = load_mmcif(cif("2rh1"), "2rh1").chains["A"]
    assert ch.uniprot_accessions() == ["P00720", "P07550"]
    assert any(d.details == "expression tag" for d in ch.seq_difs)
    assert any(d.details == "engineered mutation" for d in ch.seq_difs)
    ch.use_sifts([RefSegment("UNP", "P07550", 1, 10, 1, 10)])
    assert ch.ref_source == "sifts" and ch.uniprot_accessions() == ["P07550"]
    ch.use_sifts(None)
    assert ch.ref_source == "sifts"


def test_load_structure_rejects_non_polymer(tmp_path):
    f = tmp_path / "x.cif"
    f.write_text("data_x\n_entry.id x\n")
    with pytest.raises(StructureError):
        load_mmcif(f, "x")


def test_pdb_format_fallback_matches_mmcif(tmp_path):
    st = gemmi.read_structure(str(cif("1cns")))
    st.setup_entities()
    pdb = tmp_path / "1cns.pdb"
    st.write_pdb(str(pdb))
    m = load_structure(pdb, "1cns", fmt="pdb")
    ref = load_mmcif(cif("1cns"), "1cns").chains["A"]
    ch = m.chains["A"]
    assert m.info.fmt == "pdb" and ch.sequence == ref.sequence
    assert ch.observed.sum() == ref.observed.sum()
    assert ch.res(10).auth == ref.res(10).auth
    assert load_pdb(pdb, "1cns").chains["A"].refs == []


# ---------------------------------------------------------------- breaks
def test_breaks_detect_gaps_and_missing_atoms():
    ch = make_chain("AAAAAAAAAA", shift_after=5)
    brk = ch.breaks()
    assert brk.tolist() == [False] * 4 + [True] + [False] * 4
    ch = make_chain("AAAAA")
    ch.coords[2, 2] = np.nan  # missing C: fall back to CA-CA
    assert not ch.breaks().any()
    ch.coords[2, :] = np.nan  # residue observed but no backbone atoms
    assert ch.breaks()[1] and ch.breaks()[2]
    assert make_chain("A").breaks().size == 0


# ---------------------------------------------------------------- DSSP
def _three(c):
    return "H" if c in "HI" else ("E" if c == "E" else "-")


@pytest.mark.parametrize("pid", ["1cns", "2dch", "5nje"])
def test_dssp_matches_mkdssp_reference(pid):
    ref = json.loads((DATA / "mkdssp_reference.json").read_text())[pid]
    m = load_mmcif(cif(pid), pid)
    n = agree = 0
    for chain_id, states in ref.items():
        ch = m.chains[chain_id]
        ss = dssp.assign(ch.coords, np.array([r.mon_id == "PRO" for r in ch.residues]))
        for r in ch.residues:
            if r.observed and r.auth in states:
                n += 1
                agree += _three(ss[r.seq_id - 1]) == _three(states[r.auth])
    assert n > 100 and agree / n >= 0.99


def test_dssp_small_and_unobserved():
    assert dssp.assign(backbone(2)) == "--"
    c = backbone(6)
    c[3] = np.nan
    ss = dssp.assign(c)
    assert ss[3] == "?" and len(ss) == 6 and set(ss) <= set("-?")
    assert dssp.is_regular("H") and dssp.is_regular("E") and not dssp.is_regular("G")
    assert dssp.links(backbone(1)).size == 0


def test_hbond_matrix_blocks_are_consistent():
    m = load_mmcif(cif("5nje"), "5nje")
    ch = next(iter(m.chains.values()))
    a = dssp.hbond_matrix(ch.coords, block=512)
    b = dssp.hbond_matrix(ch.coords, block=37)
    assert (a == b).all() and a.sum() > 0
