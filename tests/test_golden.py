import pandas as pd

from linkerdb.golden import coarse, load_golden, score_golden, summary


def test_coarse_mapping():
    assert coarse("rigid_helical") == "rigid" and coarse("strand_spanning") == "rigid"
    assert coarse("flexible_disordered") == "flexible" and coarse("intermediate") == "intermediate"
    assert coarse("weird") == "undetermined"


def test_bundled_golden_set_loads():
    g = load_golden()
    assert len(g) >= 10
    assert set(g.columns) >= {"protein", "uniprot", "expected_rigidity", "confidence"}
    assert set(g["expected_rigidity"]) <= {"rigid", "flexible"}


def _canon(rows):
    return pd.DataFrame(rows)


def test_score_golden_matches_and_misses(tmp_path):
    golden = pd.DataFrame([
        {"protein": "Rigid prot", "uniprot": "P1", "expected_rigidity": "rigid", "confidence": "high", "note": ""},
        {"protein": "Flex prot", "uniprot": "P2", "expected_rigidity": "flexible", "confidence": "high", "note": ""},
        {"protein": "Mid prot", "uniprot": "P3", "expected_rigidity": "flexible", "confidence": "low", "note": ""},
        {"protein": "Missing prot", "uniprot": "P9", "expected_rigidity": "rigid", "confidence": "low", "note": ""},
    ])
    canonical = _canon([
        {"uniprot": "P1", "canonical_id": "P1:1-9", "rigidity": "rigid", "n_structures": 10,
         "rigidity_confidence": "high", "unp_lo": 1, "unp_hi": 9, "dom_dist_cv": 0.01, "elbow_angle_std": 1.0,
         "elbow_angle_range": 2.0},
        {"uniprot": "P2", "canonical_id": "P2:1-9", "rigidity": "flexible", "n_structures": 8,
         "rigidity_confidence": "high", "unp_lo": 1, "unp_hi": 9, "dom_dist_cv": 0.3, "elbow_angle_std": 20.0,
         "elbow_angle_range": 60.0},
        {"uniprot": "P3", "canonical_id": "P3:1-9", "rigidity": "intermediate", "n_structures": 5,
         "rigidity_confidence": "high", "unp_lo": 1, "unp_hi": 9, "dom_dist_cv": 0.1, "elbow_angle_std": 8.0,
         "elbow_angle_range": 15.0},
    ])
    scored = score_golden(golden, canonical)
    by = {r.protein: r for r in scored.itertuples()}
    assert by["Rigid prot"].match is True and by["Flex prot"].match is True
    assert by["Mid prot"].match is None            # intermediate = no credit either way
    assert by["Missing prot"].found is False
    text = summary(scored)
    assert "3 found" in text and "2/2 agree" in text and "NOT FOUND" in text
