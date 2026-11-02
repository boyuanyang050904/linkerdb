import json

import pytest

from linkerdb.domains import (AuthRes, CathProvider, ChainedProvider, PfamProvider, ProviderUnavailable,
                              ScopeProvider, build_provider, parse_cdf_line, parse_pdbe_pfam, parse_scope_region)
from linkerdb.domains.cath import CathFormatError, cath_paths, fetch_cath, parse_seqres_ranges
from linkerdb.domains.pfam import _resolve_overlaps
from linkerdb.download import DownloadError, NotFoundError


# ---------------------------------------------------------------- CATH CDF parsing (README examples)
def test_cdf_two_continuous_domains():
    name, doms, frags = parse_cdf_line("1chmA  D02 F00  1  A    2 - A  156 -  1  A  157 - A  402 -")
    assert name == "1chmA" and len(doms) == 2 and frags == []
    assert doms[0][0].start_auth == AuthRes("A", 2) and doms[1][0].end_auth == AuthRes("A", 402)


def test_cdf_discontinuous_domain():
    _, doms, _ = parse_cdf_line("1cnsA  D02 F00  2  A    1 - A   87 -  A  146 - A  243 -  1  A   88 - A  145 -")
    assert [(s.start_auth.num, s.end_auth.num) for s in doms[0]] == [(1, 87), (146, 243)]
    assert [(s.start_auth.num, s.end_auth.num) for s in doms[1]] == [(88, 145)]


def test_cdf_fragments_insertion_codes_and_chain_zero():
    _, doms, frags = parse_cdf_line(
        "2p070 D02 F02  2  0   16 - 0  115 -  0  231 - 0  242 -  1  0  116 - 0  230 -  "
        "0   15 A 0   15 B (2)  0  243 - 0  244 - (2)")
    assert len(doms) == 2 and len(frags) == 2
    assert frags[0].start_auth == AuthRes("0", 15, "A") and frags[0].end_auth == AuthRes("0", 15, "B")
    _, doms, _ = parse_cdf_line("1abcA D02 F00  1  A   -5 - A   10 -  1  A   11 - A   20 -")
    assert doms[0][0].start_auth.num == -5


@pytest.mark.parametrize("line", ["garbage", "1abcA D02 F00  1  A 1 - A", "1abcA D01 F00  1  A 1 - A 5 - extra",
                                  "1abcA Dxx F00"])
def test_cdf_malformed(line):
    with pytest.raises(CathFormatError):
        parse_cdf_line(line)


def test_seqres_ranges():
    assert [(s.start, s.end) for s in parse_seqres_ranges("2-78,187-208")] == [(2, 78), (187, 208)]
    assert [(s.start, s.end) for s in parse_seqres_ranges("-3-10")] == [(-3, 10)]


def test_cath_provider_from_fixture(data_dir):
    p = cath_paths(data_dir / "raw")
    prov = CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])
    assert "1cnsA" in prov.chain_keys() and "1cns" in prov.entries()
    cd = prov.get("1CNS", "A")
    assert cd.source == "CATH v4.4.0" and cd.n_domains == 2
    assert cd.domains[0].classification == "1.10.530.10" and cd.domains[0].length == 185
    assert [(s.start, s.end) for s in cd.seqres_domains[0].segments] == [(1, 87), (146, 243)]
    # unclassified domains keep their length but have no code (6r61 is not yet classified)
    d = prov.get("6r61", "A").domains[0]
    assert d.classification is None and d.length == 180
    assert prov.get("1cns", "Z") is None


def test_cath_list_overflowing_columns(tmp_path):
    cdf = tmp_path / "cdf"
    cdf.write_text("7ppeA D02 F00  1  A    1 - A  100 -  1  A  101 - A  200 -\n")
    lst = tmp_path / "list"
    lst.write_text("# comment\n7ppeA01     6     1     1     1     1     1     1     1110732   215 1.860\n"
                   "7ppeA02     3    40    50   140    99     2     1     1     8   bad 3.000\n")
    prov = CathProvider(cdf, None, lst)
    cd = prov.get("7ppe", "A")
    assert cd.domains[0].length == 215 and cd.domains[0].classification == "6.1.1.1"
    assert cd.domains[1].length is None and cd.seqres_domains is None


def test_single_domain_chains_are_ignored(tmp_path):
    cdf = tmp_path / "cdf"
    cdf.write_text("#c\n\n101mA D01 F00  1  A    0 - A  153 -\n")
    assert CathProvider(cdf).chain_keys() == []


class FakeDownloader:
    def __init__(self, fail=(), payload=None):
        self.fail = set(fail)
        self.payload = payload or {}
        self.calls = []

    def fetch(self, url, dest, validator=None):
        self.calls.append(url)
        if any(f in url for f in self.fail):
            raise DownloadError(f"blocked {url}")
        assert dest.exists(), f"test fixture missing: {dest}"
        return dest

    def get_json(self, url, method="GET", payload=None):
        self.calls.append(url)
        if url in self.payload:
            if isinstance(self.payload[url], Exception):
                raise self.payload[url]
            return self.payload[url]
        raise NotFoundError(url)


def test_fetch_cath_uses_cache(data_dir):
    fd = FakeDownloader()
    paths = fetch_cath(data_dir / "raw", fd)
    assert set(paths) == {"cdf", "seqres", "list", "unclassified"} and len(fd.calls) == 4


# ---------------------------------------------------------------- SCOPe
def test_parse_scope_region():
    out = parse_scope_region("A:1-87,A:146-243")
    assert [(c, s.start_auth.num, s.end_auth.num) for c, s in out] == [("A", 1, 87), ("A", 146, 243)]
    assert parse_scope_region("A:")[0] == ("A", None)
    (c, s), = parse_scope_region("B:-3A-120B")
    assert s.start_auth == AuthRes("B", -3, "A") and s.end_auth == AuthRes("B", 120, "B")
    assert parse_scope_region("-") == [("", None)]
    with pytest.raises(ValueError):
        parse_scope_region("A:1to5")


def test_scope_provider(data_dir, tmp_path):
    prov = ScopeProvider(data_dir / "raw/scope/dir.cla.scope.2.08-stable.txt")
    assert prov.source == "SCOPe 2.08"
    cd = prov.get("2rh1", "A")
    assert cd.n_domains == 2 and cd.domains[0].classification.startswith("f.13")
    assert prov.get("1cns", "A") is None  # SCOPe keeps chitinase as one domain
    assert "5njeA" in prov.chain_keys() and "5nje" in prov.entries()
    f = tmp_path / "cla.txt"
    f.write_text("# header\nd1abca1 1abc A:1-50 a.1.1.1 1 x\nd1abca2 1abc A:51-90 b.1.1.1 2 x\n"
                 "d1abc.1 1abc A:,B: c.1.1.1 3 x\nd1abcb1 1abc B:1x-5 c.1.1.1 4 x\nshort line\n")
    prov = ScopeProvider(f)
    assert prov.chain_keys() == ["1abcA"]


# ---------------------------------------------------------------- Pfam
def test_pfam_overlap_resolution():
    # nested region splits its container (insertion); partial overlap is cut at the midpoint
    doms = _resolve_overlaps([(57, 461, "PF00001", "7tm"), (245, 387, "PF00959", "lyz")])
    segs = {d.classification: [(s.start, s.end) for s in d.segments] for d in doms}
    assert segs == {"PF00001": [(57, 244), (388, 461)], "PF00959": [(245, 387)]}
    doms = _resolve_overlaps([(1, 100, "A", ""), (90, 200, "B", ""), (1, 100, "A", "")])
    segs = {d.classification: [(s.start, s.end) for s in d.segments] for d in doms}
    assert segs == {"A": [(1, 95)], "B": [(96, 200)]}


def test_parse_pdbe_pfam_and_provider(tmp_path):
    payload = {"1cuk": {"Pfam": {
        "PF01330": {"name": "RuvA_N", "mappings": [{"chain_id": "A", "start": {"residue_number": 1},
                                                     "end": {"residue_number": 62}}]},
        "PF14520": {"name": "HHH_5", "mappings": [{"chain_id": "A", "start": {"residue_number": 72},
                                                    "end": {"residue_number": 130}}]}}}}
    assert len(parse_pdbe_pfam(payload, "1CUK")["A"]) == 2
    url = "https://www.ebi.ac.uk/pdbe/api/mappings/pfam/1cuk"
    fd = FakeDownloader(payload={url: payload})
    prov = PfamProvider(tmp_path / "pfam", fd)
    cd = prov.get("1cuk", "A")
    assert cd.n_domains == 2 and cd.domains[0].segments[0].start == 1
    assert prov.get("1cuk", "A") is not None and fd.calls.count(url) == 1  # cached on disk
    assert prov.get("2xyz", "A") is None  # 404 -> no domains
    assert json.loads((tmp_path / "pfam/2xyz.json").read_text()) == {}


# ---------------------------------------------------------------- source selection
def test_build_provider_prefers_cath(data_dir):
    prov, notes = build_provider(data_dir / "raw", FakeDownloader())
    assert isinstance(prov, CathProvider) and notes == ["using CATH v4.4.0"]


def test_build_provider_falls_back_to_scope_then_pfam(data_dir):
    prov, notes = build_provider(data_dir / "raw", FakeDownloader(fail=["cathdb"]))
    assert isinstance(prov, ScopeProvider) and notes[0].startswith("cath unavailable")
    prov, notes = build_provider(data_dir / "raw", FakeDownloader(fail=["cathdb", "scop.berkeley"]))
    assert isinstance(prov, PfamProvider) and len(notes) == 3
    with pytest.raises(ProviderUnavailable):
        build_provider(data_dir / "raw", FakeDownloader(fail=["cathdb"]), order=("cath", "bogus"))


def test_build_provider_rejects_empty_cath(tmp_path):
    raw = tmp_path / "raw"
    (raw / "cath").mkdir(parents=True)
    for name in ("cath-domain-boundaries.txt", "cath-domain-boundaries-seqreschopping.txt",
                 "cath-domain-list.txt", "cath-unclassified-list.txt"):
        (raw / "cath" / name).write_text("")
    with pytest.raises(ProviderUnavailable):
        build_provider(raw, FakeDownloader(), order=("cath",))


def test_chained_provider(data_dir):
    p = cath_paths(data_dir / "raw")
    cath = CathProvider(p["cdf"], p["seqres"], p["list"])
    scope = ScopeProvider(data_dir / "raw/scope/dir.cla.scope.2.08-stable.txt")
    chained = ChainedProvider([scope, cath])
    assert chained.get("2rh1", "A").source.startswith("SCOPe")
    assert chained.get("1cns", "A").source.startswith("CATH")
    assert chained.get("9zzz", "A") is None and "SCOPe" in chained.source
