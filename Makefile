SHELL := /bin/zsh
UV := UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python MPLCONFIGDIR=.mplconfig uv run

.PHONY: env-check derive-program scrna cosmx cosmx-extract cosmx-analysis crc-atlas figures audit

env-check:
	$(UV) python scripts/00_environment_check.py

derive-program:
	Rscript scripts/04_reactive_egc_limma.R \
		data/raw/GSE231802/GSE231802_Counts_BulkEGCData.txt.gz \
		data/raw/GSE231802/GSE231802_series_matrix.txt.gz \
		results/tables/GSE231802_limma

scrna:
	$(UV) python scripts/01_marker_feasibility.py --dataset GSE132465 \
		--matrix data/raw/GSE132465/raw_UMI_count_matrix.txt.gz \
		--annotation data/raw/GSE132465/cell_annotation.txt.gz \
		--output results/tables/GSE132465_marker_summary.csv
	$(UV) python scripts/01_marker_feasibility.py --dataset GSE144735 \
		--matrix data/raw/GSE144735/raw_UMI_count_matrix.txt.gz \
		--annotation data/raw/GSE144735/cell_annotation.txt.gz \
		--output results/tables/GSE144735_marker_summary.csv
	$(UV) python scripts/11_scrna_threshold_sensitivity.py

cosmx-extract:
	Rscript scripts/06_cosmx_extract_schwann.R \
		data/raw/GSE303070/Processed \
		data/processed/GSE303070_schwann_targets.csv.gz
	Rscript scripts/14_cosmx_strengthening_extract.R \
		data/raw/GSE303070/Processed \
		data/processed/GSE303070_schwann_all_gene_pseudobulk.csv.gz \
		data/processed/GSE303070_schwann_identity_targets.csv.gz

cosmx-analysis:
	$(UV) python scripts/09_cosmx_schwann_pseudobulk.py
	$(UV) python scripts/16_commsbio_spatial_gradient.py

cosmx: cosmx-extract cosmx-analysis

crc-atlas:
	$(UV) python scripts/17_crc_atlas_targeted_extract.py
	$(UV) python scripts/19_crc_atlas_schwann_expression.py
	$(UV) python scripts/20_crc_atlas_external_validation.py
	$(UV) python scripts/23_crc_atlas_expression_matched_null.py

figures:
	$(UV) python scripts/12_manuscript_figures.py
	$(UV) python scripts/21_study_selection_and_source_index.py
	$(UV) python scripts/23_crc_atlas_expression_matched_null.py

audit:
	$(UV) python scripts/99_release_audit.py
