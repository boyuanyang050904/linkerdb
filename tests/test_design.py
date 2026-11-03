import numpy as np
import pandas as pd

from linkerdb import design


def _row(**kw):
    base = {"rigidity_final": "rigid", "rigidity": "rigid", "rigidity_source": "superposition",
            "linker_class": "coil", "n_structures": 20}
    base.update(kw)
    return base


def test_rigidity_trust_tiers():
    assert design.rigidity_trust(_row(linker_class="helix_spanning")) == "designed_fold"
    assert design.rigidity_trust(_row(linker_class="strand_spanning")) == "designed_fold"
    assert design.rigidity_trust(_row(n_structures=20)) == "evidence_strong"
    assert design.rigidity_trust(_row(n_structures=4)) == "evidence_weak"      # low n = possibly under-sampled
    assert design.rigidity_trust(_row(rigidity_source="proxy")) == "proxy_only"
    assert design.rigidity_trust(_row(rigidity_final="flexible")) is None      # not rigid


def test_reach_model_matches_constants():
    from linkerdb.artificial import EXTENDED_COIL, HELIX_RISE
    assert design.reach(10, "coil") == EXTENDED_COIL * 10
    assert design.reach(10, "helix") == HELIX_RISE * 10


def test_spans_thresholds():
    assert design.spans(10, "coil", 30.0) is True      # 3.5*10 = 35 >= 30
    assert design.spans(5, "coil", 30.0) is False       # 3.5*5 = 17.5 < 30
    assert design.spans(20, "helix", 25.0) is True      # 1.5*20 = 30 >= 25


def test_diversify_is_conservative_and_liability_safe():
    from linkerdb import developability
    variants = design.diversify("LGLYGKD", k=3)
    assert variants and all(len(v) == len("LGLYGKD") for v in variants)
    assert "LGLYGKD" not in variants
    s0 = developability.liability_score(developability.scan("LGLYGKD"))
    for v in variants:                                   # never worsen developability
        assert developability.liability_score(developability.scan(v)) <= s0


def test_ensemble_rigidity_rigid_vs_flexible():
    rng = np.random.RandomState(1)
    ndom = {i: rng.randn(3) * 5 for i in range(8)}
    cdom = {i: rng.randn(3) * 5 + np.array([20.0, 0, 0]) for i in range(8)}
    ref = ("ref", ndom, cdom)
    same = ("same", dict(ndom), dict(cdom))
    assert design.ensemble_rigidity([ref, same])["rigidity"] == "rigid"
    moved = ("moved", dict(ndom), {i: v + np.array([12.0, 0, 0]) for i, v in cdom.items()})
    r = design.ensemble_rigidity([ref, moved])
    assert r["rigidity"] == "flexible" and r["max_c_rmsd"] > 3.0


def _canon():
    # one of each: strong-evidence coil, weak (low-n) coil, proxy-only coil, helix fold-motif, a flexible
    return pd.DataFrame([
        {"canonical_id": "P1:1-8", "sequence": "LGLYGKD", "coil_sequence": "LGLYGKD", "linker_class": "coil",
         "rigidity_final": "rigid", "rigidity": "rigid", "rigidity_source": "superposition", "n_structures": 40,
         "length_median": 7, "is_human": True, "n_glyc": 0, "n_cys": 0, "liability_score": 0,
         "dom_dist_median": 24.0, "best_resolution": 1.8, "liability_tags": []},
        {"canonical_id": "P2:1-6", "sequence": "ARDCTT", "coil_sequence": "ARDCTT", "linker_class": "coil",
         "rigidity_final": "rigid", "rigidity": "rigid", "rigidity_source": "superposition", "n_structures": 3,
         "length_median": 6, "is_human": False, "n_glyc": 0, "n_cys": 1, "liability_score": 1,
         "dom_dist_median": 22.0, "best_resolution": 2.0, "liability_tags": ["cysteine"]},
        {"canonical_id": "P3:1-5", "sequence": "QTHVP", "coil_sequence": "QTHVP", "linker_class": "coil",
         "rigidity_final": "rigid", "rigidity": "rigid", "rigidity_source": "proxy", "n_structures": 1,
         "length_median": 5, "is_human": False, "n_glyc": 0, "n_cys": 0, "liability_score": 0,
         "dom_dist_median": 20.0, "best_resolution": 2.2, "liability_tags": []},
        {"canonical_id": "P4:1-2", "sequence": "SLDSREAFRLLS", "coil_sequence": "", "linker_class": "helix_spanning",
         "rigidity_final": "rigid", "rigidity": "rigid", "rigidity_source": "superposition", "n_structures": 6,
         "length_median": 0, "is_human": True, "n_glyc": 0, "n_cys": 0, "liability_score": 0,
         "dom_dist_median": 23.0, "best_resolution": 1.5, "liability_tags": []},
        {"canonical_id": "P5:1-9", "sequence": "GGGSGGGSG", "coil_sequence": "GGGSGGGSG", "linker_class": "coil",
         "rigidity_final": "flexible", "rigidity": "flexible", "rigidity_source": "superposition", "n_structures": 8,
         "length_median": 9, "is_human": True, "n_glyc": 0, "n_cys": 0, "liability_score": 0,
         "dom_dist_median": 24.0, "best_resolution": 1.9, "liability_tags": []},
    ])


def test_propose_filters_and_tiers():
    canon = _canon()
    df = design.propose(canon, distance=23, distance_tol=5, min_trust="evidence_weak", top=20)
    ids = list(df["candidate_id"])
    assert "P5:1-8" not in ids and all("P5" not in i for i in ids)   # flexible excluded
    assert all("P3" not in i for i in ids)                            # proxy-only below evidence_weak cutoff
    # designed_fold (continuous helix) is the top trust tier, so it ranks first;
    # the evidence_strong coil is the best splice-in option
    assert df.iloc[0]["candidate_id"] == "P4:1-2" and df.iloc[0]["trust"] == "designed_fold"
    splice = df[df["use"] == "splice"]
    assert splice.iloc[0]["candidate_id"] == "P1:1-8" and splice.iloc[0]["trust"] == "evidence_strong"
    by = {r.candidate_id: r for r in df.itertuples()}
    assert by["P1:1-8"].use == "splice" and by["P4:1-2"].use == "fold-motif"
    assert pd.isna(by["P4:1-2"].reach_A)                              # spanning: reach not meaningful


def test_propose_proxy_only_included_when_allowed():
    df = design.propose(_canon(), min_trust="proxy_only", top=20)
    assert any("P3" in i for i in df["candidate_id"])
