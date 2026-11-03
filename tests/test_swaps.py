import pandas as pd

from linkerdb.artificial import EXTENDED_COIL, HELIX_RISE
from linkerdb.swaps import load_swaps, predict, score_swaps, summary


def test_bundled_swap_set_loads():
    s = load_swaps()
    assert len(s) >= 8
    assert set(s.columns) >= {"system", "reference", "axis", "observable", "expected_kind", "expected"}
    assert set(s["axis"]) <= {"reach", "rigidity", "activity", "expression"}
    # references are flagged as unverified recollection, never presented as checked citations
    assert all("recollection" in r for r in s["reference"])


def _row(**kw):
    base = {"system": "x", "axis": "reach", "expected_kind": "", "observable": "", "expected": "",
            "la_name": "", "la_type": "", "la_len": "", "lb_name": "", "lb_type": "",
            "span_target_A": "", "tol": "", "confidence": "", "reference": "recollection"}
    base.update(kw)
    return base


def test_min_len_for_span_uses_coil_constant():
    # model lower bound = ceil(target / EXTENDED_COIL); agree if within tol of observed
    target, obs = 35, 12
    hit = predict(_row(expected_kind="min_len_for_span", expected=obs, span_target_A=target, tol=4))
    assert hit[1] is True
    assert str(int(-(-target // EXTENDED_COIL))) in hit[0]  # predicted length appears in the text
    # too tight a tolerance turns the same row into a miss
    miss = predict(_row(expected_kind="min_len_for_span", expected=obs, span_target_A=target, tol=1))
    assert miss[1] is False


def test_span_sufficiency_from_reach_bound():
    # (G4S)1 = 5 res, coil upper bound 17.5 A < 35 A -> cannot span; observed oligomer -> agree
    assert predict(_row(expected_kind="span_insufficient", la_type="coil", la_len=5, span_target_A=35))[1] is True
    # (G4S)3 = 15 res, 52.5 A >= 35 A -> can span; observed monomer -> agree
    assert predict(_row(expected_kind="span_sufficient", la_type="coil", la_len=15, span_target_A=35))[1] is True
    # a helix uses the smaller per-residue rise, so 15 res only reaches 22.5 A < 35 A
    assert HELIX_RISE * 15 < 35
    assert predict(_row(expected_kind="span_sufficient", la_type="helix", la_len=15, span_target_A=35))[1] is False


def test_length_distance_coupling_direction():
    assert predict(_row(expected_kind="length_distance_coupling", la_type="helix", expected="controls"))[1] is True
    assert predict(_row(expected_kind="length_distance_coupling", la_type="coil", expected="weak"))[1] is True
    # library would be wrong if a coil were claimed to control distance
    assert predict(_row(expected_kind="length_distance_coupling", la_type="coil", expected="controls"))[1] is False


def test_more_rigid_and_tie():
    assert predict(_row(expected_kind="more_rigid", la_type="helix", lb_type="coil", expected="a"))[1] is True
    assert predict(_row(expected_kind="more_rigid", la_type="coil", lb_type="helix", expected="a"))[1] is False
    # two coils: library cannot distinguish -> no credit either way
    assert predict(_row(expected_kind="more_rigid", la_type="coil", lb_type="coil", expected="a"))[1] is None
    # library has no call for an "other" type -> no credit
    assert predict(_row(expected_kind="more_rigid", la_type="other", lb_type="coil", expected="a"))[1] is None


def test_out_of_scope_not_scored():
    assert predict(_row(expected_kind="no_geometry_prediction", axis="expression"))[1] is None


def test_score_and_summary_split_scored_from_coverage():
    swaps = pd.DataFrame([
        _row(system="reach-hit", expected_kind="span_insufficient", la_type="coil", la_len=5, span_target_A=35),
        _row(system="rigidity-miss", axis="rigidity", expected_kind="more_rigid",
             la_type="coil", lb_type="helix", expected="a"),
        _row(system="coverage", axis="expression", expected_kind="no_geometry_prediction"),
    ])
    scored = score_swaps(swaps)
    by = {r.system: r for r in scored.itertuples()}
    assert by["reach-hit"].agree is True and by["reach-hit"].scored
    assert by["rigidity-miss"].agree is False and by["rigidity-miss"].scored
    assert by["coverage"].agree is None and not by["coverage"].scored
    text = summary(scored)
    assert "1/2" in text              # one of two scored rows agrees
    assert "SEED" in text and "out of scope" in text


def test_bundled_corpus_scores_cleanly():
    scored = score_swaps(load_swaps())
    decided = scored[scored["scored"]]
    assert len(decided) >= 5
    # the seed is hand-picked to cases the models should handle; this guards the harness,
    # not the models' accuracy (see the honesty note in swaps.py / design.md section 16)
    assert bool(decided["agree"].all())
    assert (~scored["scored"]).sum() >= 3   # expression/activity/retrieval rows are coverage-only
