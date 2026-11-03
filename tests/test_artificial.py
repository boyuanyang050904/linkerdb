from linkerdb.artificial import HELIX_RISE, artificial_table


def test_artificial_table():
    t = artificial_table()
    assert {"(G4S)1", "(EAAAK)1"} <= set(t["name"])
    # all designed to be inert -> zero liabilities
    assert (t["liability_score"] == 0).all()
    assert t.loc[t["name"].str.startswith("(G4S)"), "liability_tags"].map(len).sum() == 0
    # EAAAK rigid helix has a modelled axial reach = 1.5 A/residue; G4S flexible has none
    e = t[t["name"] == "(EAAAK)2"].iloc[0]
    assert e["class"] == "rigid_helical" and e["axial_reach_A"] == round(HELIX_RISE * 10, 1)
    import pandas as pd
    assert pd.isna(t[t["name"] == "(G4S)2"].iloc[0]["axial_reach_A"])
