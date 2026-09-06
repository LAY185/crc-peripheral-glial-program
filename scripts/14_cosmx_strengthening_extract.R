#!/usr/bin/env Rscript

# One-pass strengthening extraction for the Communications Biology gate.
# Produces all-gene Schwann pseudobulk counts plus a compact cell-level table
# for lineage and cross-lineage contamination sensitivity analyses.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "data/raw/GSE303070/Processed"
pb_output <- if (length(args) >= 2) args[[2]] else "data/processed/GSE303070_schwann_all_gene_pseudobulk.csv.gz"
identity_output <- if (length(args) >= 3) args[[3]] else "data/processed/GSE303070_schwann_identity_targets.csv.gz"

reactive_core <- c(
  "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
  "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1"
)
lineage <- c("S100B", "SOX10", "PLP1", "SLC1A3", "SOX2", "MPZ", "S100A1", "GFAP")
myeloid <- c("PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "CD68", "CD163", "APOE")
epithelial <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20")
fibroblast <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1")
wanted <- unique(c(reactive_core, lineage, myeloid, epithelial, fibroblast))

files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No RDS files found in ", input_dir)

sample_acc <- list()
identity_rows <- vector("list", length(files))

add_named <- function(a, b) {
  nm <- union(names(a), names(b))
  out <- setNames(numeric(length(nm)), nm)
  if (length(a)) out[names(a)] <- a
  if (length(b)) out[names(b)] <- out[names(b)] + b
  out
}

for (i in seq_along(files)) {
  obj <- readRDS(files[[i]])
  obj_attrs <- attributes(obj)
  md <- obj_attrs$meta.data
  assay <- attributes(obj_attrs$assays[[1]])
  layer <- attributes(assay$layers$counts)
  features <- attributes(assay$features)$dimnames[[1]]
  cells <- attributes(assay$cells)$dimnames[[1]]
  stopifnot(identical(cells, rownames(md)))

  keep_cells <- which(md$Manual_toplevel_pred == "Schwann")
  if (!length(keep_cells)) next
  mat <- new(
    "dgCMatrix", i = layer$i, p = layer$p, x = layer$x,
    Dim = layer$Dim, Dimnames = list(features, cells), factors = list()
  )

  sample_id <- as.character(md$Sample[keep_cells][1])
  donor <- as.character(md$Donor[keep_cells][1])
  sample_type <- as.character(md$Sample_type[keep_cells][1])
  mmr <- if ("MMR status" %in% names(md)) as.character(md[["MMR status"]][keep_cells][1]) else NA_character_
  lib <- md$nCount_RNA[keep_cells]
  if (all(is.na(lib))) lib <- Matrix::colSums(mat)[keep_cells]
  gene_counts <- Matrix::rowSums(mat[, keep_cells, drop = FALSE])
  names(gene_counts) <- features
  present <- setNames(rep(1, length(features)), features)

  acc <- sample_acc[[sample_id]]
  if (is.null(acc)) {
    acc <- list(
      Donor = donor, Sample_type = sample_type, MMR_status = mmr,
      n_schwann = 0, total_counts = 0, n_fovs = 0,
      counts = numeric(), fovs_present = numeric()
    )
  }
  acc$n_schwann <- acc$n_schwann + length(keep_cells)
  acc$total_counts <- acc$total_counts + sum(lib, na.rm = TRUE)
  acc$n_fovs <- acc$n_fovs + 1
  acc$counts <- add_named(acc$counts, gene_counts)
  acc$fovs_present <- add_named(acc$fovs_present, present)
  sample_acc[[sample_id]] <- acc

  present_targets <- intersect(wanted, features)
  target_counts <- t(as.matrix(mat[present_targets, keep_cells, drop = FALSE]))
  target_norm <- log1p(sweep(target_counts, 1, pmax(lib, 1), "/") * 1000)
  colnames(target_norm) <- paste0("expr_", colnames(target_norm))
  keep_md <- intersect(
    c("Barcode", "Sample", "Donor", "Sample_type", "Slide", "fov", "nCount_RNA",
      "nFeature_RNA", "niches", "Quality", "MMR status"), names(md)
  )
  ident <- cbind(md[keep_cells, keep_md, drop = FALSE], as.data.frame(target_norm))
  ident$source_file <- basename(files[[i]])
  ident$panel_n_genes <- length(features)
  identity_rows[[i]] <- ident

  if (i %% 25 == 0) message("processed ", i, "/", length(files))
}

pb_rows <- lapply(names(sample_acc), function(sample_id) {
  x <- sample_acc[[sample_id]]
  data.frame(
    Sample = sample_id,
    Donor = x$Donor,
    Sample_type = x$Sample_type,
    MMR_status = x$MMR_status,
    n_schwann = x$n_schwann,
    total_counts = x$total_counts,
    n_fovs = x$n_fovs,
    gene = names(x$counts),
    count = as.numeric(x$counts),
    fovs_present = as.numeric(x$fovs_present[names(x$counts)]),
    stringsAsFactors = FALSE
  )
})
pb <- do.call(rbind, pb_rows)

identity_names <- unique(unlist(lapply(identity_rows, names)))
identity_rows <- lapply(identity_rows[!vapply(identity_rows, is.null, logical(1))], function(x) {
  missing <- setdiff(identity_names, names(x))
  for (nm in missing) x[[nm]] <- NA
  x[, identity_names, drop = FALSE]
})
identity <- do.call(rbind, identity_rows)

dir.create(dirname(pb_output), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(pb_output, open = "wt")
write.csv(pb, con, row.names = FALSE, na = "")
close(con)
con <- gzfile(identity_output, open = "wt")
write.csv(identity, con, row.names = FALSE, na = "")
close(con)

message("wrote ", nrow(pb), " sample-gene pseudobulk rows to ", pb_output)
message("wrote ", nrow(identity), " Schwann-cell identity rows to ", identity_output)

