#!/usr/bin/env Rscript

# Extract a compact Schwann-cell expression table from the 352 GSE303070
# Seurat v5 RDS files.  We deliberately read the serialized sparse-matrix
# slots directly: this makes the extraction reproducible without installing
# the full Seurat stack merely to access public counts.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "data/raw/GSE303070/Processed"
output_file <- if (length(args) >= 2) args[[2]] else "data/processed/GSE303070_schwann_targets.csv.gz"

reactive_core <- c(
  "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
  "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1"
)
reactive_extended <- c(
  "CXCL2", "CXCL1", "STEAP4", "CCL7", "CEBPD", "NOS2", "GBP2",
  "VNN1", "SAA1", "SAA2", "GBP5", "HP", "CP", "ARG2", "IL13RA2", "FTH1"
)
context_genes <- c(
  "S100B", "GFAP", "SOX10", "PLP1", "SLC1A3", "IL1R1", "IL1B", "SPP1",
  "CD68", "APOE", "CD163", "CXCL9", "CXCL12", "TGFB1", "TNF"
)
wanted <- unique(c(reactive_core, reactive_extended, context_genes))

extract_one <- function(path) {
  obj <- readRDS(path)
  obj_attrs <- attributes(obj)
  md <- obj_attrs$meta.data
  assay <- attributes(obj_attrs$assays[[1]])
  layer <- attributes(assay$layers$counts)
  features <- attributes(assay$features)$dimnames[[1]]
  cells <- attributes(assay$cells)$dimnames[[1]]
  stopifnot(identical(cells, rownames(md)))

  keep_cells <- which(md$Manual_toplevel_pred == "Schwann")
  if (!length(keep_cells)) return(NULL)
  present <- intersect(wanted, features)
  keep_features <- match(present, features)

  mat <- new(
    "dgCMatrix", i = layer$i, p = layer$p, x = layer$x,
    Dim = layer$Dim, Dimnames = list(features, cells), factors = list()
  )
  counts <- as.matrix(mat[keep_features, keep_cells, drop = FALSE])
  counts <- t(counts)

  lib <- md$nCount_RNA[keep_cells]
  if (all(is.na(lib))) lib <- Matrix::colSums(mat)[keep_cells]
  norm <- log1p(sweep(counts, 1, pmax(lib, 1), "/") * 1000)
  colnames(norm) <- paste0("expr_", colnames(norm))

  keep_md <- c(
    "Barcode", "Sample", "Donor", "Sample_type", "Slide", "fov",
    "MMR status", "Sidedness", "Primary tumour location", "nCount_RNA",
    "nFeature_RNA", "CenterX_global_px", "CenterY_global_px",
    "CenterX_local_px", "CenterY_local_px", "niches", "Quality"
  )
  keep_md <- intersect(keep_md, names(md))
  out <- cbind(md[keep_cells, keep_md, drop = FALSE], as.data.frame(norm))
  out$source_file <- basename(path)
  out$panel_n_genes <- length(features)
  out$core_genes_present <- paste(intersect(reactive_core, features), collapse = ";")
  out
}

files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No RDS files found in ", input_dir)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  rows[[i]] <- extract_one(files[[i]])
  if (i %% 50 == 0) message("processed ", i, "/", length(files))
}

# Different panels have different columns; fill absent genes with NA rather
# than zero so that an unmeasured transcript is never treated as non-expression.
all_names <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(x) {
  missing <- setdiff(all_names, names(x))
  for (nm in missing) x[[nm]] <- NA
  x[, all_names, drop = FALSE]
})
result <- do.call(rbind, rows)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(output_file, open = "wt")
write.csv(result, con, row.names = FALSE, na = "")
close(con)

message("wrote ", nrow(result), " Schwann cells from ", length(files), " FOVs to ", output_file)
message("samples: ", length(unique(result$Sample)), "; donors: ", length(unique(result$Donor)))
