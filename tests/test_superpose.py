import numpy as np

from linkerdb.superpose import classify, mobility, relative_shift


def _rot(deg, axis="z"):
    t = np.radians(deg)
    c, s = np.cos(t), np.sin(t)
    return {"z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]),
            "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])}[axis]


def _domains(seed=0):
    rng = np.random.default_rng(seed)
    n = {i: p for i, p in enumerate(rng.normal(size=(8, 3)) * 8)}
    c = {i: p + np.array([25, 0, 0]) for i, p in enumerate(rng.normal(size=(8, 3)) * 6)}
    return n, c


def test_relative_shift_identity_and_rotation():
    n, c = _domains()
    assert relative_shift(n, c, n, c) == (0.0, 0.0)
    # rotate the C-domain 30 deg about the hinge: N unchanged -> hinge ~30, c_rmsd > 0
    c_rot = {i: _rot(30) @ p for i, p in c.items()}
    rmsd, hinge = relative_shift(n, c, n, c_rot)
    assert hinge == __import__("pytest").approx(30, abs=1) and rmsd > 1


def test_relative_shift_needs_core():
    n, c = _domains()
    small = {0: n[0], 1: n[1]}            # only 2 shared N residues (< MIN_CORE)
    assert relative_shift(n, c, small, c) is None


def test_mobility_and_classify():
    n, c = _domains()
    c30 = {i: _rot(30) @ p for i, p in c.items()}
    g = mobility([("ref", n, c), ("a", n, c), ("b", n, c30)])
    assert g["n_pairs"] == 2 and g["max_hinge"] > 25
    assert classify(g["max_c_rmsd"], g["max_hinge"]) == "flexible"
    assert mobility([("only", n, c)])["max_c_rmsd"] is None      # single structure
    assert classify(0.5, 2.0) == "rigid" and classify(None, None) == "undetermined"
    assert classify(2.0, 15.0) == "intermediate"
