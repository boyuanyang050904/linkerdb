import json

import pandas as pd
import pytest

from linkerdb.domains import CathProvider
from linkerdb.domains.cath import cath_paths
from linkerdb.download import DownloadError
from linkerdb.selection import (align_pair, attach_clusters, candidate_chains, chains_of_entries, dedup_linkers,
                                identity, map_position, rank_candidates, sample_entries)
from linkerdb.sifts import fetch_sifts, load_sifts
from linkerdb.sources import (fetch_entity_info, fetch_indices, index_paths, parse_clusters, parse_entry_type,
                              parse_resolu)


@pytest.fixture
def cath(data_dir):
    p = cath_paths(data_dir / "raw")
    return CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])


# ---------------------------------------------------------------- sources
def test_parse_indices(data_dir):
    p = index_paths(data_dir / "raw")
    res = parse_resolu(p["resolu"])
    assert res["1cns"] == pytest.approx(1.91) and res["2j4n"] is None and "idcode" not in res
    et = parse_entry_type(p["entry_type"])
    assert et["2j4n"] == ("prot", "NMR")
    cl = parse_clusters(p["clusters"])
    assert "1CNS_1" in cl
    assert set(parse_clusters(p["clusters"], {"1CNS_1"})) == {"1CNS_1"}


class FakeDL:
    def __init__(self, entries=None, fail_batches=()):
        self.entries = entries or {}
        self.fail = set(fail_batches)
        self.posts = 0
        self.fetched = []

    def fetch(self, url, dest, validator=None):
        self.fetched.append(url)
        return dest

    def get_json(self, url, method="GET", payload=None):
        self.posts += 1
        if self.posts in self.fail:
            raise DownloadError("boom")
        ids = payload["variables"]["ids"]
        return {"data": {"entries": [self.entries[i] for i in ids if i in self.entries] + [None]}}


def _gql(pid, chains, unp):
    return {"rcsb_id": pid, "rcsb_entry_info": {"resolution_combined": [2.1, 1.9], "deposited_model_count": 1},
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "polymer_entities": [{"rcsb_id": f"{pid}_1", "entity_poly": {"type": "polypeptide(L)"},
                                  "rcsb_polymer_entity_container_identifiers": {
                                      "entity_id": "1", "auth_asym_ids": chains, "uniprot_ids": unp}}]}


def test_fetch_entity_info_batches_cache_and_resume(tmp_path):
    cache = tmp_path / "e.jsonl"
    dl = FakeDL({"1ABC": _gql("1ABC", ["A"], ["P1"]), "2ABC": _gql("2ABC", ["B"], [])}, fail_batches={2})
    info = fetch_entity_info(["1abc", "2abc", "3abc"], cache, dl, batch=1)
    assert info["1abc"]["resolution"] == 1.9 and info["1abc"]["entities"][0]["auth_asym_ids"] == ["A"]
    assert "2abc" not in info  # its batch failed
    assert info["3abc"]["missing"] is True  # obsolete entry: remembered as missing
    with open(cache, "a") as fh:
        fh.write('{"pdb_id": "trunc')  # interrupted write
    dl2 = FakeDL({"2ABC": _gql("2ABC", ["B"], [])})
    info = fetch_entity_info(["1abc", "2abc", "3abc"], cache, dl2, batch=10)
    assert dl2.posts == 1 and set(info) == {"1abc", "2abc", "3abc"}


def test_fetch_indices_and_sifts(data_dir):
    dl = FakeDL()
    paths = fetch_indices(data_dir / "raw", dl)
    assert len(dl.fetched) == 3 and paths["clusters"].name == "clusters-by-entity-40.txt"
    sifts = load_sifts(fetch_sifts(data_dir / "raw", dl))
    assert sifts["2rh1"]["A"][1].accession == "P00720"
    assert set(load_sifts(fetch_sifts(data_dir / "raw", dl), {"1cns"})) <= {"1cns"}


# ---------------------------------------------------------------- selection
def test_sample_entries_is_deterministic(cath, data_dir):
    et = parse_entry_type(index_paths(data_dir / "raw")["entry_type"])
    a = sample_entries(cath, et, n_diffraction=3, n_nmr=1, seed=1)
    assert a == sample_entries(cath, et, n_diffraction=3, n_nmr=1, seed=1)
    assert len(a) == 4 and "2j4n" in a
    chains = chains_of_entries(cath, a)
    assert all(k in a for k in chains) and all(chains.values())


def test_candidates_clusters_and_ranking(cath, data_dir):
    p = index_paths(data_dir / "raw")
    cands = candidate_chains(cath, parse_resolu(p["resolu"]), parse_entry_type(p["entry_type"]), 2.5)
    assert set(cands["pdb_id"]) == {"1cns", "2rh1", "6r61", "2dch", "5nje"}  # 2f7z is 3.0 Å, 2j4n NMR
    assert set(candidate_chains(cath, parse_resolu(p["resolu"]), parse_entry_type(p["entry_type"]), 1.9)["pdb_id"]) \
        == {"6r61"}
    ents = {}
    for line in open(data_dir / "raw/rcsb/entities.jsonl"):
        rec = json.loads(line)
        ents[rec["pdb_id"]] = rec
    clusters = parse_clusters(p["clusters"])
    n_dom = {k: cath.get(k[:4], k[4:]).n_domains for k in cands["pdb_id"] + cands["chain"]}
    df = attach_clusters(cands, ents, clusters, n_dom)
    assert df.loc[df.pdb_id == "2rh1", "n_uniprot"].item() == 2
    ranked = rank_candidates(df)
    # 1cns A and B are copies of one entity -> one chain kept
    assert (ranked["pdb_id"] == "1cns").sum() == 1
    assert ranked.groupby("cluster")["rank"].min().eq(0).all()


def test_attach_clusters_unmapped_and_chain_zero():
    cands = pd.DataFrame([{"pdb_id": "1abc", "chain": "0", "resolution": 2.0, "method_class": "diffraction"},
                          {"pdb_id": "2abc", "chain": "A", "resolution": 2.0, "method_class": "diffraction"},
                          {"pdb_id": "3abc", "chain": "A", "resolution": 2.0, "method_class": "diffraction"}])
    ents = {"1abc": {"entities": [{"entity": "1ABC_1", "type": "polypeptide(L)", "auth_asym_ids": ["A"],
                                   "uniprot_ids": ["P1"]}]},
            "2abc": {"entities": [{"entity": "2ABC_1", "type": "polypeptide(L)", "auth_asym_ids": ["A"],
                                   "uniprot_ids": []}]}}
    df = attach_clusters(cands, ents, {"1ABC_1": 7}, {})
    assert df["cluster"].tolist() == ["c7", "single:2ABC_1", "nomap:3abcA"]


def test_alignment_helpers():
    a = "MKTAYIAKQRQISFVKSHFSRQ"
    assert identity(a, a) == 1.0 and identity(a, "") == 0.0
    ident, aln = align_pair("XXXX" + a, a)
    assert ident == 1.0
    assert map_position(aln, 5) == 1 and map_position(aln, 1) == 1  # overhang maps to nearest block
    assert identity("MKTAUYIAKQ", "MKTACYIAKQ") == 0.9  # selenocysteine (U) is not in BLOSUM62: no crash


def _lk(lid, pdb, chain, acc, lo, hi, cath_n="1.1.1.1", cath_c="2.2.2.2", res=2.0, ss=(10, 12)):
    return {"linker_id": lid, "pdb_id": pdb, "auth_asym_id": chain, "uniprot": acc, "unp_lo": lo, "unp_hi": hi,
            "cath_n": cath_n, "cath_c": cath_c, "resolution": res, "junction_index": 1,
            "ss_start": ss[0], "ss_end": ss[1]}


def test_dedup_linkers():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLKLEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQCGVWPAACRESVKEALV"
    mutant = seq[:100] + "W" * 20 + seq[120:]
    unrelated = "".join(reversed(seq))
    df = pd.DataFrame([
        _lk("a_A_1", "aaaa", "A", "P1", 100, 110, res=1.5),
        _lk("a_A_2", "aaaa", "A", "P1", 105, 115, res=1.5),          # same chain: never compared
        _lk("b_A_1", "bbbb", "A", "P1", 108, 120, res=2.0),          # same UniProt region as a_A_1
        _lk("c_A_1", "cccc", "A", "P9", 1, 5, res=2.2, ss=(10, 12)),   # homolog of a (different accession)
        _lk("d_A_1", "dddd", "A", "P8", 1, 5, res=2.3, ss=(200, 202)),  # homolog chain but other junction
        _lk("e_A_1", "eeee", "A", "P7", 1, 5, res=2.4),               # unrelated sequence
        _lk("f_A_1", "ffff", "A", None, None, None, cath_n=None, res=2.5),
    ])
    seqs = {("aaaa", "A"): seq, ("bbbb", "A"): seq, ("cccc", "A"): mutant, ("dddd", "A"): mutant,
            ("eeee", "A"): unrelated, ("ffff", "A"): seq}
    out = dedup_linkers(df, ["resolution", "linker_id"], 0.4, chain_seqs=seqs)
    dropped = out.attrs["dropped"]
    assert set(out["linker_id"]) == {"a_A_1", "a_A_2", "d_A_1", "e_A_1", "f_A_1"}
    assert dropped["b_A_1"].startswith("same_uniprot_region:a_A_1")
    assert dropped["c_A_1"].startswith("homologous_chain:a_A_1")
    # without chain sequences only the UniProt rule applies
    assert set(dedup_linkers(df, ["resolution"])["linker_id"]) == set(df["linker_id"]) - {"b_A_1"}
