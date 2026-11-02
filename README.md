# linkerdb

Systematic extraction of inter-domain **linkers** of natural multi-domain proteins from the PDB.
This repository currently contains the **data acquisition layer**: structure download, domain-boundary
parsing, linker extraction and dataset construction. Design, rules and known limitations are in
[`docs/design.md`](docs/design.md).

## Install

```bash
pip install -e ".[test]"
```

## Usage

```bash
linkerdb fetch-meta                          # CATH v4.4.0, wwPDB indices, RCSB clusters (cached in ./data)
linkerdb sample --out results/sample50       # 40 diffraction + 10 NMR random entries, summary printed
linkerdb build  --out results/full           # <= 2.5 Å, 40 % non-redundant, writes parquet files
linkerdb summary results/full/linkers.parquet
linkerdb show 5t6i_A_1 13pk_A_1              # readable view of single linkers (sequence / CATH / DSSP tracks)
```

`--data-dir` (or `LINKERDB_DATA`) sets the cache location; `--workers` sets the process count.
Every download is cached and resumable (HTTP Range); every processed chain is check-pointed under
`data/work/<run>/entries/`, so an interrupted `build` continues where it stopped.

## Outputs

| file | content |
|---|---|
| `linkers.parquet` | final non-redundant, QC-passed linkers (one row per domain junction) |
| `linkers_all.parquet` | every junction of every processed chain, with `qc_flags`, `qc_pass`, `dedup_dropped` |
| `chains.parquet` | per-chain status (`ok`, `excluded_chimera`, `numbering_disagree`, ...) and checks |
| `candidates.parquet` | candidate chains with their 40 % cluster and rank (build only) |
| `run_info.json` | data versions, rule parameters and counts of the run |

Key columns: `gap_*` = residues CATH leaves unassigned between two domains (strict definition, often empty);
`ss_*` = structural linker (coil stretch between the last SSE of the N-terminal and the first SSE of the
C-terminal domain, DSSP-based). Positions are SEQRES indices (`label_seq_id`); `*_auth` are author numbers.

## Tests

```bash
pytest --cov=linkerdb          # offline; fixtures in tests/data (regenerate: scripts/make_test_fixtures.py)
```
