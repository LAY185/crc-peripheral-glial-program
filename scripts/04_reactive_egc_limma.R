args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 04_reactive_egc_limma.R counts.gz series_matrix.gz output_dir")
}

counts_path <- args[[1]]
series_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(limma))

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
series_lines <- readLines(gzfile(series_path))
title_line <- series_lines[grepl("^!Sample_title", series_lines)][1]
titles <- strsplit(title_line, "\t")[[1]][-1]
titles <- gsub('^"|"$', "", titles)
description_line <- series_lines[grepl("^!Sample_description", series_lines)][1]
descriptions <- strsplit(description_line, "\t")[[1]][-1]
descriptions <- gsub('^"|"$', "", descriptions)
if (length(titles) != ncol(counts)) stop("Sample title/count mismatch")

extract <- regexec("EGCs, (.*), ([0-9]+)h, rep([0-9]+)", titles)
parts <- regmatches(titles, extract)
meta_by_geo_order <- data.frame(
  sample = descriptions,
  treatment = vapply(parts, `[[`, character(1), 2),
  time = vapply(parts, `[[`, character(1), 3),
  replicate = as.integer(vapply(parts, `[[`, character(1), 4)),
  stringsAsFactors = FALSE
)
if (!setequal(meta_by_geo_order$sample, colnames(counts))) {
  stop("Sample descriptions do not map to count-matrix columns")
}
meta <- meta_by_geo_order[match(colnames(counts), meta_by_geo_order$sample), ]
meta$treatment <- gsub(" conditioned medium", "_CM", meta$treatment)
meta$group <- factor(paste0("h", meta$time, "_", gsub(" ", "_", meta$treatment)))
write.csv(meta, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

keep <- rowSums(counts >= 10) >= 4
filtered <- counts[keep, , drop = FALSE]
libsize <- colSums(filtered)
log_cpm <- log2(t(t(filtered) / libsize) * 1e6 + 0.5)
design <- model.matrix(~0 + meta$group)
colnames(design) <- levels(meta$group)
fit <- lmFit(log_cpm, design)

contrast_text <- c(
  tumor_vs_healthy_6h = "h6_tumor_CM-h6_healthy_CM",
  tumor_vs_healthy_12h = "h12_tumor_CM-h12_healthy_CM",
  tumor_vs_healthy_24h = "h24_tumor_CM-h24_healthy_CM",
  tumor_vs_unstim_6h = "h6_tumor_CM-h6_unstimulated",
  tumor_vs_unstim_12h = "h12_tumor_CM-h12_unstimulated",
  tumor_vs_unstim_24h = "h24_tumor_CM-h24_unstimulated"
)
contrasts <- makeContrasts(contrasts = contrast_text, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE)

tables <- list()
for (nm in colnames(contrasts)) {
  tab <- topTable(fit2, coef = nm, number = Inf, sort.by = "P")
  tab$gene <- rownames(tab)
  tab$contrast <- nm
  tables[[nm]] <- tab
  write.csv(tab, file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
}

healthy_names <- grep("^tumor_vs_healthy", names(tables), value = TRUE)
merged <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = TRUE),
  lapply(healthy_names, function(nm) {
    z <- tables[[nm]][, c("gene", "logFC", "adj.P.Val")]
    colnames(z)[-1] <- paste0(c("logFC_", "FDR_"), sub("tumor_vs_healthy_", "", nm))
    z
  })
)
lfc_cols <- grep("^logFC_", colnames(merged), value = TRUE)
fdr_cols <- grep("^FDR_", colnames(merged), value = TRUE)
merged$n_up_fdr05 <- rowSums(merged[, lfc_cols, drop = FALSE] > 0.5 & merged[, fdr_cols, drop = FALSE] < 0.05)
merged$n_down_fdr05 <- rowSums(merged[, lfc_cols, drop = FALSE] < -0.5 & merged[, fdr_cols, drop = FALSE] < 0.05)
merged$mean_logFC <- rowMeans(merged[, lfc_cols, drop = FALSE], na.rm = TRUE)
merged <- merged[order(-merged$n_up_fdr05, -merged$mean_logFC), ]
write.csv(merged, file.path(output_dir, "reactive_egc_consensus.csv"), row.names = FALSE)

strict <- merged[merged$n_up_fdr05 >= 2, ]
writeLines(strict$gene, file.path(output_dir, "reactive_egc_strict_genes.txt"))

cat("Retained genes:", nrow(filtered), "\n")
cat("Strict reactive-EGC genes (up in >=2 time points):", nrow(strict), "\n")
cat("Top candidates:\n")
print(head(merged[, c("gene", lfc_cols, fdr_cols, "n_up_fdr05", "mean_logFC")], 25))
