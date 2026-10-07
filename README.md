# How differently generated evaluation matrices shape biomedical model performance, model comparison, drug ranking and candidate selection: a cross-ecosystem analysis of LINCS, DepMap and GTEx — reproduction materials

This directory has two parts:

- `github/` — the **code and documentation** pushed to GitHub (version control; easy to file issues/PRs)
- `zenodo/` — the **derived data, figures, manuscript and checksums** uploaded to Zenodo (immutable archive; gets a DOI)

## Shortest reproduction path (reviewer's view)

1. Download `zenodo/source_data/` from the Zenodo records: frozen prediction matrices, DepMap-derived data,
   per-sample GTEx scores, and figure source data.
2. Clone the code repository from GitHub and run `pip install -r requirements.txt`.
3. Run `run_all.sh`: the scripts recompute Tables 1–5, Figures 1–7 and Supplementary Tables S1–S26b in order,
   and compare each key number against `EXPECTED_OUTPUTS.json` (tolerance 5×10⁻⁴).
4. For the parts that must be rerun from the raw public data (LINCS level-5, DepMap 21Q2, GTEx v11), see
   `DATA_SOURCES.md`. For the LINCS side, the exported reference matrices and the `inst_id` lists are provided,
   so reviewers can either use the matrices directly or rebuild them deterministically from the public products.

## What we cannot redistribute (documented in DATA_SOURCES.md)

- LINCS L1000 raw / level-5 products (GEO GSE92742, SigCom LINCS, DCIC 2021): public but large — download and rebuild paths given.
- DepMap 21Q2 expression / dependency matrices: public (CC-BY 4.0) — download and subset-build paths given.
- GTEx v11 raw expression: public (dbGaP / GTEx portal) — download path given.
- MSigDB C2:CGP three disease signatures: restricted by the MSigDB license — this package gives only the set names and version;
  download the gene lists yourself under that license.
