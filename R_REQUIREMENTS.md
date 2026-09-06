# R environment

The analyses used R 4.5.1 with these CRAN/Bioconductor packages:

- `limma`
- `Matrix`

Install the packages before running the R targets:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
install.packages("Matrix")
```
