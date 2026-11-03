import numpy as np
import pandas as pd
import pytest

from linkerdb import cli, developability
from linkerdb.aggregate import _merge_intervals, _rigidity, aggregate
from linkerdb.linkers import chain_ss, extract_linkers, map_domains
from linkerdb.retrieve import format_hits, query_canonical
from linkerdb.structure import load_mmcif

from .conftest import DATA, make_chain
from .test_linkers import mapped, two_domains


# ---------------------------------------------------------------- developability
def test_developability_scan():
    c = developability.scan("NASNGTDGSMWCRАRKR".replace("А", ""))  # strip a stray char
    assert c["deamid_ng"] == 1 and c["isomer_dg"] == 1 and c["n_cys"] == 1
    assert developability.scan("NPS")["n_glyc"] == 0            # X != P
    assert developability.scan("NISNIT")["n_glyc"] == 2         # overlapping sequons both counted
    assert developability.scan("RSRR")["furin"] == 1
    tags = developability.summary_flags(developability.scan("NITDG"))
    assert "N-glycosylation" in tags and "isomerization(DG)" in tags
    assert developability.liability_score(developability.scan("NIT")) > developability.liability_score(
        developability.scan("AAAA"))
    assert developability.scan("") == {k: 0 for k in developability.LIABILITY_KEYS}


# ---------------------------------------------------------------- geometry fields
def test_geometry_fields_present():
    ch, m = two_domains(n=60, cut=30)
    (rec,) = extract_linkers(ch, m, "H" * 28 + "-" * 4 + "H" * 28)
    assert rec["dom_dist"] is not None and rec["elbow_angle"] is not None
    # extended straight chain: the two domain centroids and the junction are collinear-ish
    assert rec["elbow_angle"] > 150


# ---------------------------------------------------------------- aggregation
def test_merge_intervals():
    assert _merge_intervals([(1, 10, 0), (5, 15, 1), (30, 40, 2)]) == [[0, 1], [2]]
    assert _merge_intervals([(1, 5, 0), (6, 9, 1)]) == [[0], [1]]  # touching but not overlapping


def _inst(pdb, acc, lo, hi, dom_dist, ss="coil", miss=0, ss_len=8, e2e=20.0, ang=120.0):
    return dict(linker_id=f"{pdb}_A_1", pdb_id=pdb, auth_asym_id="A", uniprot=acc, unp_lo=lo, unp_hi=hi,
                ss_length=ss_len, ss_sequence="GSGSGSGS"[:ss_len] or "G", junction_ss=ss, junction_type="sequential",
                end_to_end=e2e, dom_dist=dom_dist, elbow_angle=ang, n_missing=miss, resolution=2.0,
                cath_n="3.1.1.1", cath_c="2.1.1.1", taxid="9606", qc_pass=True, qc_flags=[], method="X-RAY",
                description="test")


def test_rigidity_uses_boundary_independent_distance():
    # same junction, but linker boundary drawn differently across structures:
    # end-to-end swings wildly while the domains stay put -> must read as rigid, not flexible
    rows = [_inst(p, "P1", 100, 110, dom_dist=30.0 + d, ss_len=L, e2e=E)
            for p, d, L, E in [("1aaa", 0.1, 9, 25.0), ("1bbb", -0.1, 3, 7.0), ("1ccc", 0.2, 9, 25.5),
                               ("1ddd", 0.0, 4, 9.0)]]
    out = aggregate(pd.DataFrame(rows))
    assert len(out) == 1
    r = out.iloc[0]
    assert r["rigidity"] == "rigid" and r["rigidity_confidence"] == "high"
    assert r["dom_dist_cv"] < 0.05 and r["end_to_end_cv"] > 0.3   # e2e noise is reported but not used


def test_rigidity_flexible_and_helix():
    flex = pd.DataFrame([_inst(p, "P2", 50, 70, dom_dist=dd) for p, dd in
                         [("2aaa", 20.0), ("2bbb", 35.0), ("2ccc", 55.0)]])
    assert aggregate(flex).iloc[0]["rigidity"] == "flexible"
    helix = pd.DataFrame([_inst(p, "P3", 10, 30, dom_dist=40.0, ss="helix_spanning", ss_len=0)
                          for p in ("3a", "3b", "3c")])
    # ss_len 0 rows are dropped as non-samples; give them length
    helix = pd.DataFrame([_inst(p, "P3", 10, 30, dom_dist=40.0 + i, ss="helix_spanning")
                          for i, p in enumerate(("3a", "3b", "3c"))])
    assert aggregate(helix).iloc[0]["rigidity"] == "rigid_helical"


def test_aggregate_groups_and_singletons():
    rows = [_inst("1aaa", "P1", 100, 110, 30.0), _inst("1bbb", "P1", 102, 112, 30.1),  # one group
            _inst("1ccc", "P1", 300, 320, 40.0)]                                         # separate region
    ungrouped = _inst("2xxx", None, None, None, 25.0)
    ungrouped["unp_lo"] = None
    out = aggregate(pd.DataFrame(rows + [ungrouped]))
    assert len(out) == 3
    g = out[out["grouped"]]
    assert set(g["n_structures"]) == {2, 1} and (~out["grouped"]).sum() == 1
    assert out["is_human"].sum() == 3  # taxid 9606 for the P1 rows
    assert out.loc[~out["grouped"], "canonical_id"].iloc[0] == "2xxx_A_1"


def test_aggregate_keeps_spanning_drops_failed():
    # a zero-length (SSE-spanning) junction is kept and labelled; a qc-failed one is dropped
    rows = pd.DataFrame([dict(_inst("1aaa", "P1", 100, 110, 30.0, ss_len=0, ss="helix_spanning"),
                              ss_sequence="", flank_n_seq="AAAEEEHHHHHH", flank_c_seq="HHHHHHKKKAAA"),
                         dict(_inst("1bbb", "P1", 100, 110, 30.0), qc_pass=False)])
    out = aggregate(rows)
    assert len(out) == 1
    r = out.iloc[0]
    assert r["linker_class"] == "helix_spanning" and r["length_median"] == 0
    assert r["coil_sequence"] == "" and r["sequence"]        # window from the flanking helix
    assert r["rigidity"] == "rigid_helical"


# ---------------------------------------------------------------- retrieval
@pytest.fixture
def canon():
    rows = []
    for i in range(6):
        rows.append({"canonical_id": f"P{i}:1-10", "uniprot": f"P{i}", "description": "d", "length_median": 6 + i,
                     "rigidity": "rigid" if i % 2 else "flexible", "rigidity_confidence": "high",
                     "dom_dist_median": 20.0 + 3 * i, "end_to_end_median": 18.0, "end_to_end_cv": 0.1,
                     "elbow_angle_median": 120.0, "n_structures": 10 - i, "is_human": i < 3,
                     "n_glyc": i % 2, "n_cys": 0, "liability_score": i, "liability_tags": [],
                     "best_resolution": 1.5, "sequence": "GSGSGS"})
    return pd.DataFrame(rows)


def test_query_filters_and_ranks(canon):
    hits = query_canonical(canon, distance=26.0, distance_tol=3.5)
    assert (abs(hits["dom_dist_median"] - 26.0) <= 3.5).all()
    assert hits.iloc[0]["canonical_id"] == "P2:1-10"       # closest to 26 (dom_dist 26)
    assert query_canonical(canon, human_only=True)["is_human"].all()
    assert (query_canonical(canon, rigidity="rigid")["rigidity"] == "rigid").all()
    assert (query_canonical(canon, no_glyc=True)["n_glyc"] == 0).all()
    assert len(query_canonical(canon, top=2)) == 2
    assert "no linkers match" in format_hits(query_canonical(canon, distance=999, distance_tol=1))
    assert "P2:1-10" in format_hits(query_canonical(canon, distance=26.0, distance_tol=3.5))


# ---------------------------------------------------------------- end to end, offline
def test_build_and_query_offline(data_dir, tmp_path, capsys):
    common = ["--data-dir", str(data_dir), "--workers", "1"]
    out = tmp_path / "db"
    assert cli.main(common + ["build", "--run-name", "dbt", "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "canonical linkers:" in text
    canonical = pd.read_parquet(out / "linkers_canonical.parquet")
    instances = pd.read_parquet(out / "linkers_instances.parquet")
    assert len(canonical) > 0 and "rigidity" in canonical and "n_glyc" in canonical
    assert "cluster40" in instances and "dom_dist" in instances
    # 2rh1 is an engineered chimera: its junctions must be excluded from the canonical layer
    assert not canonical["uniprot"].isin(["P00720"]).any()
    assert cli.main(common + ["query", "--canonical", str(out / "linkers_canonical.parquet"), "--top", "5"]) == 0
    assert cli.main(common + ["query", "--canonical", str(out / "linkers_canonical.parquet"),
                              "--distance", "30", "--distance-tol", "40", "--human"]) == 0
    assert len(capsys.readouterr().out) > 0
