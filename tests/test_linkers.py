import pytest

from linkerdb.config import Rules
from linkerdb.domains import AuthRes, CathProvider, ChainDomains, Domain, Segment
from linkerdb.domains.cath import cath_paths
from linkerdb.linkers import (ChainSkip, MappedChain, MappedDomain, chain_ss, drop_minor_segments,
                              extract_linkers, junctions, map_domains, natural_problem)
from linkerdb.structure import RefSegment, SeqDif, load_mmcif

from .conftest import DATA, make_chain


def auth_dom(did, ranges, chain="A", code="1.1.1.1", length=None):
    return Domain(did, [Segment(start_auth=AuthRes(chain, a), end_auth=AuthRes(chain, b)) for a, b in ranges],
                  code, length)


def seq_dom(did, ranges, code="1.1.1.1"):
    return Domain(did, [Segment(start=a, end=b) for a, b in ranges], code)


def mapped(*doms):
    return MappedChain([MappedDomain(d, c, None, segs) for d, c, segs in doms], [], "agree")


# ---------------------------------------------------------------- mapping / numbering checks
def test_map_domains_agree_and_length_check():
    ch = make_chain("A" * 40, observed="1" * 35 + "0" * 5, auth_offset=10)
    cd = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(11, 30)], length=20), auth_dom("d2", [(31, 45)], length=9)],
                      seqres_domains=[seq_dom("d1", [(1, 20)]), seq_dom("d2", [(21, 35)])])
    m = map_domains(ch, cd)
    assert m.numbering_check == "agree" and m.domains[1].segments == [(21, 35)]
    assert "domain_length_mismatch" in m.flags  # d2: 15 observed vs 9 reported


def test_map_domains_resolves_disagreement_by_observed_endpoints():
    # SEQRES grew at the N-terminus (unmodelled tag): CATH seqres chopping is stale
    ch = make_chain("A" * 50, observed="0" * 10 + "1" * 40, auth_offset=-10)
    cd = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 40)])],
                      seqres_domains=[seq_dom("d1", [(1, 20)]), seq_dom("d2", [(21, 40)])])
    m = map_domains(ch, cd)
    assert m.numbering_check == "resolved_cdf" and m.domains[0].segments == [(11, 30)]
    assert "numbering_resolved" in m.flags
    # the reverse case: author numbering stale, SEQRES chopping right
    cd2 = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(-9, 0)]), auth_dom("d2", [(21, 30)])],
                       seqres_domains=[seq_dom("d1", [(11, 30)]), seq_dom("d2", [(31, 50)])])
    assert map_domains(ch, cd2).numbering_check == "resolved_seqres"


def test_map_domains_unresolvable_and_failed():
    ch = make_chain("A" * 40, auth_offset=0)
    cd = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 40)])],
                      seqres_domains=[seq_dom("d1", [(1, 18)]), seq_dom("d2", [(19, 40)])])
    with pytest.raises(ChainSkip) as exc:
        map_domains(ch, cd)
    assert exc.value.reason == "numbering_disagree" and "d1" in exc.value.detail
    bad = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(500, 600)])], seqres_domains=[seq_dom("d1", [(0, 3)])])
    with pytest.raises(ChainSkip) as exc:
        map_domains(ch, bad)
    assert exc.value.reason == "mapping_failed"


def test_map_domains_single_sources_and_unobserved_endpoint():
    ch = make_chain("A" * 40, observed="1" * 39 + "0")
    cdf_only = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 39)])],
                            seqres_domains=[seq_dom("d1", [(1, 20)]), seq_dom("d2", [(21, 99)])])
    assert map_domains(ch, cdf_only).numbering_check == "cdf_only"
    seq_only = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 999)])],
                            seqres_domains=[seq_dom("d1", [(1, 20)]), seq_dom("d2", [(21, 39)])])
    assert map_domains(ch, seq_only).numbering_check == "seqres_only"
    scope_like = ChainDomains("9xyz", "A", "SCOPe", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 40)])])
    assert map_domains(ch, scope_like).numbering_check == "single_source"  # endpoint check only for CATH
    pfam_like = ChainDomains("9xyz", "A", "Pfam", [seq_dom("p1", [(1, 20)], None), seq_dom("p2", [(25, 40)], None)])
    m = map_domains(ch, pfam_like)
    assert m.numbering_check == "single_source" and "cath_unclassified" in m.flags
    cath_bad = ChainDomains("9xyz", "A", "CATH", [auth_dom("d1", [(1, 20)]), auth_dom("d2", [(21, 40)])],
                            seqres_domains=[seq_dom("d1", [(1, 20)]), seq_dom("d2", [(21, 40)])])
    with pytest.raises(ChainSkip, match="unmodelled"):
        map_domains(ch, cath_bad)


# ---------------------------------------------------------------- junctions
def test_junction_types():
    m = mapped(("a", "x", [(1, 50), (101, 150)]), ("b", "y", [(51, 100)]), ("c", "z", [(151, 200)]))
    jx, flags = junctions(m)
    assert [(j.n_dom.domain_id, j.c_dom.domain_id, j.kind) for j in jx] == [
        ("a", "b", "insertion"), ("b", "a", "insertion"), ("a", "c", "sequential")]
    assert flags == []
    m = mapped(("a", "x", [(1, 50), (101, 150)]), ("b", "y", [(51, 100), (151, 200)]))
    assert {j.kind for j in junctions(m)[0]} == {"interleaved"}
    # internal gap of a single domain is not a junction; overlapping segments are flagged
    m = mapped(("a", "x", [(1, 50), (60, 100)]), ("b", "y", [(95, 150)]))
    jx, flags = junctions(m)
    assert len(jx) == 1 and jx[0].kind == "interleaved" and flags == ["overlapping_segments"]


def test_drop_minor_segments():
    m = mapped(("a", "x", [(1, 10), (21, 30)]), ("b", "y", [(11, 20)]))
    ss = "HHHH------" + "----------" + "EEEE------"
    m2, dropped = drop_minor_segments(m, ss)
    assert dropped == [("b", 11, 20)] and len(m2.domains) == 1 and "minor_segment_dropped" in m2.flags


# ---------------------------------------------------------------- linker extraction
def two_domains(n=60, cut=30, **kw):
    ch = make_chain("G" * n, **kw)
    return ch, mapped(("d1", "1.10.1.1", [(1, cut)]), ("d2", "2.20.2.2", [(cut + 1, n)]))


def test_structural_linker_extends_over_coil():
    ch, m = two_domains()
    ss = "H" * 24 + "-" * 6 + "-" * 4 + "E" * 26
    (rec,) = extract_linkers(ch, m, ss)
    assert (rec["gap_length"], rec["gap_start"], rec["gap_end"]) == (0, 31, 30)
    assert (rec["ss_start"], rec["ss_end"], rec["ss_length"]) == (25, 34, 10)
    assert (rec["ext_n"], rec["ext_c"]) == (6, 4) and rec["junction_ss"] == "coil"
    assert rec["ss_sequence"] == "G" * 10 and rec["flank_n_dssp"] == "H" * 10 and rec["flank_c_dssp"] == "E" * 10
    assert rec["qc_pass"] and rec["qc_flags"] == [] and rec["uniprot"] == "P12345"
    assert (rec["unp_start"], rec["unp_end"]) == (125, 134)
    assert rec["ca_junction"] == pytest.approx(3.8) and rec["end_to_end"] == pytest.approx(9 * 3.8)
    assert rec["ss_start_auth"] == "25" and len(rec["window_seq"]) == 10 + 2 + 40


def test_extension_cap():
    ch, m = two_domains()
    (rec,) = extract_linkers(ch, m, "H" + "-" * 58 + "H", Rules(max_ext=5))
    assert rec["ss_length"] == 10 and {"capped_n", "capped_c"} <= set(rec["qc_flags"])
    assert rec["qc_pass"]  # capping is informational


@pytest.mark.parametrize("ss,expected", [
    ("H" * 60, "helix_spanning"),
    ("E" * 60, "strand_spanning"),
    ("H" * 30 + "E" * 30, "sse_abutting"),
])
def test_zero_length_junction_classes(ss, expected):
    ch, m = two_domains()
    (rec,) = extract_linkers(ch, m, ss)
    assert rec["ss_length"] == 0 and rec["junction_ss"] == expected
    assert rec["ss_sequence"] == "" and rec["ss_start"] == rec["ss_end"] + 1 and rec["end_to_end"] is None


def test_gap_linker_and_long_gap():
    ch = make_chain("A" * 100)
    m = mapped(("d1", "x", [(1, 30)]), ("d2", "y", [(36, 100)]))
    (rec,) = extract_linkers(ch, m, "H" * 100)
    assert rec["gap_length"] == 5 and rec["gap_sequence"] == "AAAAA" and rec["ss_length"] == 5
    m = mapped(("d1", "x", [(1, 30)]), ("d2", "y", [(70, 100)]))
    (rec,) = extract_linkers(ch, m, "H" * 100)
    assert "long_gap" in rec["qc_flags"] and not rec["qc_pass"]


def test_missing_residues_and_disorder():
    obs = "1" * 28 + "0" * 6 + "1" * 26
    ch, m = two_domains(observed=obs)
    ss = "H" * 28 + "?" * 6 + "H" * 26
    (rec,) = extract_linkers(ch, m, ss)
    assert rec["ss_length"] == 6 and rec["n_missing"] == 6 and rec["frac_observed"] == 0
    assert {"fully_disordered", "missing_residues"} <= set(rec["qc_flags"]) and rec["qc_pass"]
    assert rec["ca_junction"] is None


def test_chain_break_and_geometry_violation():
    # mild break (stretched C-N bond, CA-CA still 3.8 Å): chain_break only
    ch, m = two_domains()
    ch.coords[33, 0, 1] += 1.7  # move N of residue 34 away from C of 33 (C-N ~2.2 Å)
    ss = "H" * 28 + "-" * 8 + "H" * 24
    (rec,) = extract_linkers(ch, m, ss)
    assert rec["n_breaks"] == 1 and "chain_break" in rec["qc_flags"]
    assert "geometry_violation" not in rec["qc_flags"] and rec["qc_pass"]
    # modelled neighbours 14 Å apart inside the linker: impossible for a continuous chain
    ch, m = two_domains(shift_after=33)
    (rec,) = extract_linkers(ch, m, ss)
    assert {"chain_break", "geometry_violation"} <= set(rec["qc_flags"]) and not rec["qc_pass"]
    ch, m = two_domains(shift_after=30)  # consecutive junction residues 10 Å apart
    (rec,) = extract_linkers(ch, m, "H" * 60)
    assert "geometry_violation" in rec["qc_flags"] and not rec["qc_pass"]


def test_geometry_tolerance_accepts_real_spread():
    ch, m = two_domains()
    ch.coords[30:, :, 0] += 0.08  # CA-CA 3.88 Å, as observed in real structures
    (rec,) = extract_linkers(ch, m, "H" * 60)
    assert "geometry_violation" not in rec["qc_flags"]


def test_nonstandard_residues():
    ch, m = two_domains(mon={31: "MSE", 32: "SEP"})
    (rec,) = extract_linkers(ch, m, "H" * 28 + "-" * 6 + "H" * 26)
    assert rec["n_nonstandard"] == 2 and rec["nonstandard"] == ["MSE", "SEP"]
    assert "nonstandard_residue" in rec["qc_flags"] and rec["residues_3"][2] == "MSE"


def test_uniprot_checks():
    ss = "H" * 28 + "-" * 4 + "H" * 28
    # residues 31-32 not covered by UniProt (e.g. an engineered insertion)
    refs = [RefSegment("UNP", "P1", 1, 30, 1, 30), RefSegment("UNP", "P1", 33, 60, 31, 58)]
    ch, m = two_domains(refs=refs)
    (rec,) = extract_linkers(ch, m, ss)
    assert "unp_unmapped" in rec["qc_flags"] and not rec["qc_pass"]
    # contiguous SEQRES but a jump in UniProt numbering (deleted natural segment)
    refs = [RefSegment("UNP", "P1", 1, 30, 1, 30), RefSegment("UNP", "P1", 31, 60, 51, 80)]
    ch, m = two_domains(refs=refs)
    assert "unp_discontinuous" in extract_linkers(ch, m, ss)[0]["qc_flags"]
    # two different proteins meet at the junction (fusion)
    refs = [RefSegment("UNP", "P1", 1, 30, 1, 30), RefSegment("UNP", "P2", 31, 60, 1, 30)]
    ch, m = two_domains(refs=refs)
    (rec,) = extract_linkers(ch, m, ss)
    assert "unp_discontinuous" in rec["qc_flags"] and rec["uniprot"] is None


def test_seq_dif_annotations():
    ss = "H" * 28 + "-" * 4 + "H" * 28
    ch, m = two_domains(difs=[SeqDif(30, "GLY", "ALA", "engineered mutation"), SeqDif(5, "GLY", "", "conflict")])
    (rec,) = extract_linkers(ch, m, ss)
    assert rec["n_mutations"] == 1 and rec["qc_pass"]
    ch, m = two_domains(difs=[SeqDif(31, "GLY", "", "linker")])
    assert "engineered_linker" in extract_linkers(ch, m, ss)[0]["qc_flags"]


def test_extra_flags_and_natural_problem():
    ch, m = two_domains()
    (rec,) = extract_linkers(ch, m, "H" * 60, extra_flags=["chimera"])
    assert not rec["qc_pass"]
    assert natural_problem(make_chain("AAA")) is None
    assert natural_problem(make_chain("AAA", polymer_type="polydeoxyribonucleotide"))[0] == "not_protein"
    assert natural_problem(make_chain("AAA", refs=[RefSegment("PDB", "9XYZ", 1, 3, 1, 3)])) == ("no_uniprot", "PDB")
    two = [RefSegment("UNP", "P1", 1, 2, 1, 2), RefSegment("UNP", "P2", 3, 3, 1, 1)]
    assert natural_problem(make_chain("AAA", refs=two)) == ("chimera", "P1,P2")


# ---------------------------------------------------------------- real structures
@pytest.fixture(scope="module")
def cath():
    p = cath_paths(DATA / "raw")
    return CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])


def _real(cath, pid, chain):
    ch = load_mmcif(DATA / "raw" / "mmcif" / pid[1:3] / f"{pid}.cif.gz", pid).chains[chain]
    m = map_domains(ch, cath.get(pid, chain))
    return ch, m, extract_linkers(ch, m, chain_ss(ch))


def test_real_insertion_domain_1cns(cath):
    ch, m, recs = _real(cath, "1cns", "A")
    assert m.numbering_check == "agree"
    assert [(r["domain_n"], r["domain_c"], r["junction_type"]) for r in recs] == [
        ("1cnsA01", "1cnsA02", "insertion"), ("1cnsA02", "1cnsA01", "insertion")]
    assert all(r["qc_pass"] for r in recs) and all(r["ss_start"] <= 88 <= r["ss_end"] + 1 for r in recs[:1])


def test_real_numbering_resolution_6r61(cath):
    ch, m, recs = _real(cath, "6r61", "A")
    assert m.numbering_check == "resolved_cdf" and m.domains[0].segments == [(70, 249)]
    assert recs and recs[0]["domain_n"] == "6r61A01"


def test_real_fusion_2rh1_is_flagged(cath):
    ch, m, recs = _real(cath, "2rh1", "A")
    assert natural_problem(ch)[0] == "chimera"
    assert all("unp_discontinuous" in r["qc_flags"] for r in recs)


def test_geometry_violation_across_missing_residues():
    # cleaved loop (like serpin 3nda): 2 unmodelled residues but the modelled ends are 40 Å apart
    obs = "1" * 31 + "00" + "1" * 27
    ch, m = two_domains(observed=obs, shift_after=33, cut=31)
    ch.coords[33:, :, 0] += 30.0
    ss = "H" * 28 + "----" + "??" + "--" + "H" * 24
    (rec,) = extract_linkers(ch, m, ss)
    assert rec["ca_junction"] is None  # junction residue unmodelled: the old check could not see it
    assert "geometry_violation" in rec["qc_flags"] and not rec["qc_pass"]
