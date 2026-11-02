"""Regenerate the offline test fixtures in tests/data from a populated data cache.

Usage: python scripts/make_test_fixtures.py [DATA_DIR]

Copies a handful of real mmCIF files and extracts the matching lines from the
CATH / SCOPe / SIFTS / wwPDB / RCSB metadata files, so that the whole pipeline
(including the CLI) can be tested without network access.
"""

import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

ENTRIES = ["1cns", "2rh1", "6r61", "2dch", "5nje", "2f7z", "2j4n"]
DSSP_CHECK = ["1cns", "2dch", "5nje"]  # reference mkdssp strings stored for these

root = Path(__file__).resolve().parents[1]
data = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "data"
raw = data / "raw"
out = root / "tests" / "data"
(out / "raw").mkdir(parents=True, exist_ok=True)
want = set(ENTRIES)


def copy_lines(src: Path, dst: Path, keep, header=lambda line: False, gz=False):
    dst.parent.mkdir(parents=True, exist_ok=True)
    op = gzip.open if gz else open
    with op(src, "rt") as fi, op(dst, "wt") as fo:
        for line in fi:
            if header(line) or keep(line):
                fo.write(line)


for pid in ENTRIES:
    src = raw / "mmcif" / pid[1:3] / f"{pid}.cif.gz"
    dst = out / "raw" / "mmcif" / pid[1:3] / f"{pid}.cif.gz"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)

cath = raw / "cath"
copy_lines(cath / "cath-domain-boundaries.txt", out / "raw/cath/cath-domain-boundaries.txt",
           lambda l: l[:4] in want, header=lambda l: l.startswith("#"))
for name in ("cath-domain-boundaries-seqreschopping.txt", "cath-domain-list.txt", "cath-unclassified-list.txt"):
    copy_lines(cath / name, out / "raw/cath" / name, lambda l: l[:4] in want, header=lambda l: l.startswith("#"))
copy_lines(raw / "scope/dir.cla.scope.2.08-stable.txt", out / "raw/scope/dir.cla.scope.2.08-stable.txt",
           lambda l: l.split("\t")[1:2] and l.split("\t")[1] in want, header=lambda l: l.startswith("#"))
copy_lines(raw / "sifts/pdb_chain_uniprot.tsv.gz", out / "raw/sifts/pdb_chain_uniprot.tsv.gz",
           lambda l: l.split("\t")[0] in want, header=lambda l: l.startswith("#") or l.startswith("PDB"), gz=True)
idx = raw / "index"
copy_lines(idx / "resolu.idx", out / "raw/index/resolu.idx", lambda l: l[:4].lower() in want,
           header=lambda l: not (len(l) > 4 and l[4] == "\t"))
copy_lines(idx / "pdb_entry_type.txt", out / "raw/index/pdb_entry_type.txt", lambda l: l[:4] in want)
copy_lines(idx / "clusters-by-entity-40.txt", out / "raw/index/clusters-by-entity-40.txt",
           lambda l: any(t.split("_")[0].lower() in want for t in l.split()))
ents = raw / "rcsb" / "entities.jsonl"
copy_lines(ents, out / "raw/rcsb/entities.jsonl", lambda l: json.loads(l)["pdb_id"] in want)

# reference secondary structure from mkdssp (if installed)
ref = {}
if shutil.which("mkdssp"):
    for pid in DSSP_CHECK:
        tmp = Path("/tmp") / f"{pid}.cif"
        with gzip.open(out / "raw" / "mmcif" / pid[1:3] / f"{pid}.cif.gz", "rb") as fi, open(tmp, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        dssp_out = tmp.with_suffix(".dssp")
        subprocess.run(["mkdssp", "--output-format", "dssp", str(tmp), str(dssp_out)], check=True,
                       capture_output=True)
        started = False
        for line in open(dssp_out):
            if line.startswith("  #  RESIDUE"):
                started = True
                continue
            if not started or line[13] == "!":
                continue
            ref.setdefault(pid, {}).setdefault(line[11], {})[f"{int(line[5:10])}{line[10].strip()}"] = \
                line[16] if line[16] != " " else "-"
    (out / "mkdssp_reference.json").write_text(json.dumps(ref))
print("fixtures written to", out)
