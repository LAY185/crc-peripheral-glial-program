#!/usr/bin/env python3
"""Sensitivity of scRNA reactive-core effects to the minimum glial-cell count."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/tables/cross_platform_meta/scrna_cell_threshold_sensitivity.csv"
CORE = [
    "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
    "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
]


def hedges_g(tumour: pd.Series, normal: pd.Series) -> tuple[float, float, float]:
    t = tumour.to_numpy(float)
    n = normal.to_numpy(float)
    nt, nn = len(t), len(n)
    if nt < 2 or nn < 2:
        return np.nan, np.nan, np.nan
    pooled = np.sqrt(((nt - 1) * t.var(ddof=1) + (nn - 1) * n.var(ddof=1)) / (nt + nn - 2))
    if pooled == 0:
        return np.nan, np.nan, np.nan
    g = (1 - 3 / (4 * (nt + nn) - 9)) * (t.mean() - n.mean()) / pooled
    se = np.sqrt((nt + nn) / (nt * nn) + g**2 / (2 * (nt + nn - 2)))
    return g, g - 1.96 * se, g + 1.96 * se


def main() -> None:
    rows = []
    for dataset in ["GSE132465", "GSE144735"]:
        d = pd.read_csv(ROOT / f"results/tables/{dataset}_marker_summary.csv")
        d = d[d["Cell_subtype"].eq("Enteric glial cells") & d["Class"].isin(["Tumor", "Normal"])].copy()
        genes = [f"{gene}__mean" for gene in CORE if f"{gene}__mean" in d]
        d["score"] = d[genes].mean(axis=1)
        for threshold in [1, 3, 5, 10]:
            x = d[d["n_cells"] >= threshold]
            t = x.loc[x["Class"].eq("Tumor"), "score"]
            n = x.loc[x["Class"].eq("Normal"), "score"]
            g, lo, hi = hedges_g(t, n)
            p = mannwhitneyu(t, n, alternative="two-sided").pvalue if len(t) and len(n) else np.nan
            rows.append({
                "dataset": dataset,
                "minimum_glial_cells": threshold,
                "n_tumour": len(t),
                "n_normal": len(n),
                "hedges_g": g,
                "ci_low": lo,
                "ci_high": hi,
                "mannwhitney_p": p,
            })
    result = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT, index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
