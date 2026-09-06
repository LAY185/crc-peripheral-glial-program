# CRC peripheral-glial program

This repository reproduces the analyses and figures for "Cross-platform transcriptomics identifies a recurrent tumour-associated peripheral-glial module in colorectal cancer."

The study derives a 33-gene response from tumour-conditioned mouse enteric glia, fixes a 13-gene cross-platform scaffold, and tests a `TIMP1-IFITM1-CHI3L1` module in human colorectal cancer single-cell and spatial datasets. The primary CRC Atlas analysis uses paired patients as biological replicates.

## Contents

- `scripts/`: analysis, targeted public-data extraction, statistical testing, and figure generation.
- `config/genesets.yaml`: fixed 33-gene, 13-gene, and three-gene lists.
- `config/mouse_human_scaffold_mapping.csv`: locked mouse-to-human mapping for the 13-gene scaffold.
- `results/source_data/`: panel-level source data for Fig. 1 to Fig. 4 and Supplementary Fig. S1 to Fig. S2.
- `results/tables/`: key intermediate summary tables and statistical outputs.
- `results/manuscript_figures/`: final Fig. 1 to Fig. 4 and Supplementary Fig. S1 exports.
- `results/figures/`: Supplementary Fig. S2 exports.
- `data/README.md`: public accessions, required filenames, and local layout.

## Public data

The analysis reuses public data from GEO accessions `GSE231802`, `GSE132465`, `GSE144735`, and `GSE303070`, plus the integrated CRC Atlas at `https://crc.icbi.at/`. This repository does not redistribute the downloaded raw matrices, 352 CosMx RDS files, the 32.7-GB CRC Atlas h5ad file, or API caches. See `data/README.md` for the expected input paths.

## Software environment

The locked Python environment targets Python 3.11 or 3.12. The manuscript analyses used Python 3.12 and R 4.5.1.

Install [`uv`](https://docs.astral.sh/uv/) and create the Python environment:

```bash
uv sync --frozen
uv run python scripts/00_environment_check.py
```

Install the R packages listed in `R_REQUIREMENTS.md` before running the R scripts.

## Reproduction order

Place the source files under `data/raw/` as described in `data/README.md`. The commands below write analysis-ready extracts to the ignored `data/processed/` directory and regenerate the committed summary tables and figures.

```bash
# 1. Derive the mouse perturbation program.
make derive-program

# 2. Build supportive single-cell summaries.
make scrna

# 3. Extract CosMx Schwann cells, run paired pseudobulk and spatial analyses.
make cosmx

# 4. Query the CRC Atlas, run paired validation and the matched-set null.
make crc-atlas

# 5. Rebuild source-data indices and all final figures.
make figures
```

The complete computation includes large public downloads, 352 CosMx input files, remote byte-range requests, 10,000 matched random sets, 200,000 bootstrap samples, and 2,000 spatial permutations. Run the targets separately if compute or network access is limited.

## Figure commands

```bash
# Fig. 1 to Fig. 4
uv run python scripts/12_manuscript_figures.py

# Supplementary Fig. S1 and source-data manifest
uv run python scripts/21_study_selection_and_source_index.py

# Supplementary Fig. S2 and its sensitivity tables
uv run python scripts/23_crc_atlas_expression_matched_null.py
```

The main figure script reads committed intermediate tables, so readers can regenerate the submitted figures without redownloading the large public datasets.

## Statistical design

- GSE231802 uses repository-designated replicate libraries for program derivation.
- CRC Atlas uses paired patients within studies and a hierarchical bootstrap over studies and patients.
- GSE303070 uses paired donors for pseudobulk validation and tumour specimens for the spatial-boundary analysis.
- GSE132465 and GSE144735 provide descriptive patient-tissue summaries because their paired and unpaired contributions do not support one common inferential test.
- Cells contribute to aggregates and never serve as independent biological replicates.

## Verification

Run the release audit before citing or archiving a modified version:

```bash
uv run python scripts/99_release_audit.py
```

The audit checks fixed gene-set cardinality, required source-data files, ignored raw-data paths, oversized files, and common credential or absolute-path patterns.

## Licence and citation

The analysis software is available under the MIT licence. Public source datasets retain their original repository terms. Cite the archived Zenodo release using `CITATION.cff`.
