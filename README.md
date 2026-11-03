# linkerdb

A database of inter-domain **linkers** from natural multi-domain proteins in the PDB, built for
**fusion-protein / biologics engineering** (bispecifics, scFv/Fc fusions, tandem domains): find a
natural linker that holds two domains at a needed distance with a needed rigidity, with its
developability liabilities flagged. Design, rules, the run record and known limitations are in
[`docs/design.md`](docs/design.md) (the database layer is §13).

## Install

```bash
pip install -e ".[test]"
```

## Usage

```bash
linkerdb fetch-meta                    # CATH v4.4.0, wwPDB indices, RCSB clusters, SIFTS (cached in ./data)
linkerdb build --out results/db        # all multi-domain chains <= 2.5 Å -> instances + canonical linkers
linkerdb query --distance 25 --rigidity rigid --human --no-glyc --min-structures 3
linkerdb score-golden                  # check rigidity calls against the seed golden set
linkerdb show 5t6i_A_1 13pk_A_1        # readable view of one junction (sequence / CATH / DSSP tracks)
linkerdb sample --out results/sample50 # 50-entry validation sample (for checking the extraction rules)
```

`--data-dir` (or `LINKERDB_DATA`) sets the cache; `--workers` the process count. Every download is
cached and resumable (HTTP Range); every processed chain is check-pointed under `data/work/<run>/entries/`,
so an interrupted `build` continues where it stopped.

### Querying

`linkerdb query` searches `linkers_canonical.parquet`:

```bash
linkerdb query --canonical results/db/linkers_canonical.parquet \
    --distance 25 --distance-tol 6 --rigidity rigid --human --no-glyc --min-structures 3 --top 20
```

`--distance` targets the inter-domain centroid distance (Å); `--rigidity rigid|flexible|any`;
`--min-structures N` keeps only linkers with enough structures for a high-confidence rigidity call;
`--no-glyc` / `--no-cys` / `--max-liability` drop developability-risky linkers.

## Two layers

- **`linkers_instances.parquet`** — every domain junction in every chain (one row per junction), with
  the structural-linker sequence/boundaries, per-structure geometry (`dom_dist`, `elbow_angle`), QC flags
  and the 40 % cluster id. This is the raw evidence; redundancy is kept on purpose.
- **`linkers_canonical.parquet`** — instances grouped by UniProt region into one row per natural linker,
  with consensus length/sequence, a **rigidity** call (`rigid` / `flexible` / `intermediate` /
  `rigid_helical` / …) and its confidence, geometry spread, developability liabilities, species, and how
  many structures back it. This is what `query` reads.

Rigidity is judged from how much the **domain-centroid distance** varies across structures (boundary-
independent), so it needs several structures of the same protein; `rigidity_confidence` and
`n_structures` say how much to trust it. Non-redundancy is a column (`cluster40`), not a deletion.

**Scope:** this narrows and ranks candidates and gives experimentally-observed starting conformations;
it does not predict that a linker will work in a particular construct — model and test the shortlist.

Other files: `chains.parquet` (per-chain status), `candidates.parquet` (candidate chains + clusters),
`run_info.json` (versions, rule parameters, counts).

## Status

`results/db/` is the database (CATH v4.4.0, ≤ 2.5 Å, built 2026-11-03): 114,433 junctions →
106,869 QC-passed instances → **27,266 canonical linkers** (13,419 UniProt accessions, 3,280 human);
7,914 (29 %) have ≥ 3 structures for a high-confidence rigidity call. By class: coil 24,161,
helix_spanning 1,532, strand_spanning 1,337, sse_abutting 236. `python scripts/check_database.py results/db`
re-runs the invariant checks.

`results/sample50/` is the validation sample for the extraction rules (`docs/design.md` §8–§11). The
first-stage non-redundant build and its linker-level dedup were retired in v3 (the database groups by
UniProt region instead of deleting redundant linkers). Roadmap (golden-set calibration, artificial-linker
comparison, ECOD/antibody coverage, AlphaFold layer) is `docs/design.md` §13.7.

## Tests

```bash
pytest --cov=linkerdb          # offline; fixtures in tests/data (regenerate: scripts/make_test_fixtures.py)
```
