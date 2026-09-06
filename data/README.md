# Public input data

This repository does not redistribute the large public datasets or local API caches. Download the source files from their original repositories and preserve the directory layout below.

## GEO accessions

| Accession | Role | Expected local files |
|---|---|---|
| GSE231802 | Mouse enteric-glial perturbation and 33-gene program derivation | `data/raw/GSE231802/GSE231802_Counts_BulkEGCData.txt.gz`; `data/raw/GSE231802/GSE231802_series_matrix.txt.gz` |
| GSE132465 | Supportive human single-cell analysis | `data/raw/GSE132465/raw_UMI_count_matrix.txt.gz`; `data/raw/GSE132465/cell_annotation.txt.gz` |
| GSE144735 | Supportive human single-cell analysis | `data/raw/GSE144735/raw_UMI_count_matrix.txt.gz`; `data/raw/GSE144735/cell_annotation.txt.gz` |
| GSE303070 | Paired CosMx validation and spatial-boundary analysis | 352 processed RDS files under `data/raw/GSE303070/Processed/`; author-provided all-cell metadata at `data/raw/GSE303070/All_cells_metadata.csv.gz` |

GEO landing pages follow `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=ACCESSION`.

## CRC Atlas

The scripts use the public integrated CRC Atlas without downloading its complete 32.7-GB h5ad file:

- h5ad endpoint: `https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad`
- CELLxGENE expression endpoint: `https://crc.icbi.at/atlas/api/v0.2/data/var`

The pipeline writes targeted metadata and expression extracts under `data/processed/` and caches remote responses under ignored `data/raw/` paths. It excludes the `Lee_2020_Nat_Genet` study from the primary CRC Atlas validation because GSE132465 and GSE144735 support separate analyses.

## Rights and provenance

Each source dataset retains the licence and reuse terms of its original repository and publication. The MIT licence in this repository applies to the software, not to third-party raw data.
