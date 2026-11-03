import json

import pandas as pd
import pytest

from linkerdb import cli
from linkerdb.domains import CathProvider
from linkerdb.domains.cath import cath_paths
from linkerdb.pipeline import Run, process_entry, write_outputs
from linkerdb.sifts import load_sifts, sifts_path
from linkerdb.stats import composition, show_linker, summarize
from linkerdb.structure import load_mmcif
from linkerdb.linkers import chain_ss, extract_linkers, map_domains

from .conftest import DATA


@pytest.fixture
def cath(data_dir):
    p = cath_paths(data_dir / "raw")
    return CathProvider(p["cdf"], p["seqres"], p["list"], unclassified_path=p["unclassified"])


def cif(settings, pid):
    return settings.structure_path(pid)


# ---------------------------------------------------------------- process_entry
def test_process_entry_one_per_entity(settings, cath):
    chains = [("A", cath.get("1cns", "A")), ("B", cath.get("1cns", "B"))]
    out = process_entry("1cns", chains, cif(settings, "1cns"), "cif", settings.rules, one_per_entity=True)
    status = {c["chain"]: c["status"] for c in out["chains"]}
    assert sorted(status.values()) == ["duplicate_entity_copy", "ok"]
    ok = next(c for c in out["chains"] if c["status"] == "ok")
    assert ok["numbering_check"] == "agree" and ok["n_junctions"] == 2 and ok["dssp_vs_annotation"] > 0.85
    assert len(out["linkers"]) == 2 and all(r["domain_source"] == "CATH v4.4.0" for r in out["linkers"])


def test_process_entry_failures(settings, cath, tmp_path):
    cd = cath.get("1cns", "A")
    out = process_entry("1cns", [("A", cd)], None, "cif", settings.rules, fetch_error="obsolete")
    assert out["chains"][0]["status"] == "download_failed" and out["chains"][0]["detail"] == "obsolete"
    junk = tmp_path / "junk.cif"
    junk.write_text("this is not a cif file")
    out = process_entry("1cns", [("A", cd)], junk, "cif", settings.rules)
    assert out["chains"][0]["status"] == "parse_failed"
    out = process_entry("1cns", [("Q", cd)], cif(settings, "1cns"), "cif", settings.rules)
    assert out["chains"][0]["status"] == "chain_not_found"
    # CATH chain '0' (legacy blank chain id) is remapped when the entry has one protein chain
    out = process_entry("2dch", [("0", cath.get("2dch", "X"))], cif(settings, "2dch"), "cif", settings.rules)
    assert out["chains"][0]["status"] == "ok" and "chain0_remap" in out["chains"][0]["chain_flags"]


def test_process_entry_natural_and_numbering(settings, cath):
    out = process_entry("2rh1", [("A", cath.get("2rh1", "A"))], cif(settings, "2rh1"), "cif", settings.rules)
    assert out["chains"][0]["status"] == "excluded_chimera"
    assert out["linkers"] and not any(r["qc_pass"] for r in out["linkers"])
    sifts = load_sifts(sifts_path(settings.raw_dir), {"6r61"})["6r61"]
    out = process_entry("6r61", [("A", cath.get("6r61", "A"))], cif(settings, "6r61"), "cif", settings.rules,
                        unp=sifts)
    c = out["chains"][0]
    assert c["status"] == "ok" and c["numbering_check"] == "resolved_cdf" and c["unp_source"] == "sifts"


def test_numbering_disagreement_is_reported(settings, cath):
    cd = cath.get("1cns", "A")
    cd.seqres_domains[0].segments[0].end = 80  # make the two CATH files disagree
    out = process_entry("1cns", [("A", cd)], cif(settings, "1cns"), "cif", settings.rules)
    assert out["chains"][0]["status"] == "numbering_disagree"


# ---------------------------------------------------------------- Run / checkpoints
def test_run_checkpoints_and_resume(settings, cath, monkeypatch):
    jobs = [("1cns", [("A", cath.get("1cns", "A"))]), ("5nje", [("A", cath.get("5nje", "A"))]),
            ("9zzz", [("A", cath.get("1cns", "A"))])]
    import linkerdb.pipeline as pl
    from linkerdb.download import NotFoundError

    def fake_fetch(ids, settings_, downloader=None, workers=None):
        return {i.upper(): ((settings_.structure_path(i), "cif") if i != "9zzz" else NotFoundError("gone"))
                for i in ids}

    monkeypatch.setattr(pl, "fetch_structures", fake_fetch)
    run = Run(settings, "t")
    chains, linkers = run.execute(jobs, workers=1)
    assert set(chains["status"]) == {"ok", "download_failed"} and len(linkers) >= 3
    assert run.done("1cns", "A") and (run.ckpt / "1cns_A.json").exists()

    calls = []
    monkeypatch.setattr(pl, "_worker", lambda a: calls.append(a) or pl.process_entry(*a[:7], a[7]))
    chains2, linkers2 = run.execute(jobs, workers=1)
    assert calls == []  # everything resumed from checkpoints
    assert len(chains2) == len(chains) and len(linkers2) == len(linkers)
    run.execute(jobs[:1], workers=1, force=True)
    assert len(calls) == 1


def test_write_outputs(tmp_path):
    chains = pd.DataFrame([{"pdb_id": "1abc", "chain": "A", "status": "ok"}])
    linkers = pd.DataFrame([{"linker_id": "x", "qc_flags": ["a"], "nonstandard": [], "qc_pass": True}])
    paths = write_outputs(tmp_path / "o", chains, linkers, linkers, {"n": 1})
    assert set(paths) == {"chains", "linkers_all", "linkers"}
    assert pd.read_parquet(paths["linkers"])["qc_flags"].iloc[0].tolist() == ["a"]
    assert json.loads((tmp_path / "o/run_info.json").read_text())["n"] == 1


# ---------------------------------------------------------------- stats
def test_summarize_and_show(settings, cath):
    out = process_entry("5nje", [("A", cath.get("5nje", "A"))], cif(settings, "5nje"), "cif", settings.rules)
    chains, linkers = pd.DataFrame(out["chains"]), pd.DataFrame(out["linkers"])
    text = summarize(chains, linkers, chains["sequence"], title="t")
    assert "qc_pass" in text and "composition" in text and "top CATH pairs" in text
    assert "no linkers" in summarize(chains, pd.DataFrame())
    assert composition([]) == {} and abs(sum(composition(["ACDA"]).values()) - 1) < 1e-9
    ch = load_mmcif(cif(settings, "5nje"), "5nje").chains["A"]
    ss = chain_ss(ch)
    rec = extract_linkers(ch, map_domains(ch, cath.get("5nje", "A")), ss)[0]
    view = show_linker(ch, rec, ss)
    assert rec["linker_id"] in view and "DSSP" in view and "CATH" in view


# ---------------------------------------------------------------- CLI (offline: everything is cached)
def test_cli_end_to_end(data_dir, tmp_path, capsys):
    common = ["--data-dir", str(data_dir), "--workers", "1"]
    assert cli.main(common + ["fetch-meta"]) == 0
    out = tmp_path / "sample"
    assert cli.main(common + ["sample", "--n-diffraction", "3", "--n-nmr", "1", "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "sample: 4 entries" in text and (out / "linkers.parquet").exists()

    db = tmp_path / "db"
    assert cli.main(common + ["build", "--out", str(db)]) == 0
    text = capsys.readouterr().out
    assert "canonical linkers:" in text
    instances = pd.read_parquet(db / "linkers_instances.parquet")
    chains = pd.read_parquet(db / "chains.parquet")
    info = json.loads((db / "run_info.json").read_text())
    assert info["domain_source"] == "CATH v4.4.0" and info["mode"] == "database"
    ok = instances[instances["qc_pass"]]
    assert (ok["resolution"] <= 2.5).all()
    assert "2rh1" not in set(ok["pdb_id"])  # engineered fusion never passes QC
    assert chains.loc[chains.pdb_id == "2rh1", "status"].item() == "excluded_chimera"
    assert (db / "candidates.parquet").exists()

    assert cli.main(common + ["summary", str(db / "linkers_instances.parquet")]) == 0
    assert "qc_pass" in capsys.readouterr().out
    lid = ok["linker_id"].iloc[0]
    assert cli.main(common + ["show", lid, "1cns_A_9", "1cns_Q_1", "--json"]) == 0
    text = capsys.readouterr().out
    assert lid in text and "no such junction" in text and "chain not in" in text


def test_cli_build_limit_and_parallel(data_dir, tmp_path, capsys):
    common = ["--data-dir", str(data_dir), "--workers", "2"]
    assert cli.main(common + ["build", "--limit", "2", "--max-resolution", "2.1", "--out", str(tmp_path / "f")]) == 0
    assert "canonical linkers:" in capsys.readouterr().out
