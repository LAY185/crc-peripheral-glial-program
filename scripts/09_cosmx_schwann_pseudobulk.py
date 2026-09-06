#!/usr/bin/env python3
"""Pseudobulk validation of the mouse-derived reactive-glia program in CosMx.

Counts are recovered exactly from the saved log1p(count/library*1000) values,
summed over author-annotated Schwann cells, and normalized only after summing.
Paired tumour-normal donors are the primary analysis unit.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/processed/GSE303070_schwann_targets.csv.gz"
OUT = ROOT / "results/tables/GSE303070_cosmx"
FIG = ROOT / "results/figures"
CORE = [
    "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
    "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(INPUT, low_memory=False)
    for gene in CORE:
        # Inverse of the extractor's normalization; rounding removes FP noise.
        d[f"count_{gene}"] = np.rint(
            np.expm1(d[f"expr_{gene}"]) * d["nCount_RNA"] / 1000
        ).astype(int)

    group = ["Sample", "Donor", "Sample_type", "MMR status"]
    named = {gene: (f"count_{gene}", "sum") for gene in CORE}
    pb = d.groupby(group, dropna=False).agg(
        n_schwann=("Barcode", "size"), total_counts=("nCount_RNA", "sum"), **named
    ).reset_index()
    for gene in CORE:
        pb[f"logCPM_{gene}"] = np.log1p(pb[gene] / pb["total_counts"] * 1e6)
    score_cols = [f"logCPM_{g}" for g in CORE]
    pb["reactive_core_score"] = pb[score_cols].mean(axis=1)
    pb.to_csv(OUT / "schwann_pseudobulk_scores.csv", index=False)

    wide_score = pb.pivot(index="Donor", columns="Sample_type", values="reactive_core_score").dropna()
    score_stat, score_p = wilcoxon(wide_score["T"], wide_score["N"])
    delta_score = (wide_score["T"] - wide_score["N"]).to_numpy()
    rng = np.random.default_rng(20260904)
    bootstrap_medians = np.median(
        rng.choice(delta_score, size=(200_000, len(delta_score)), replace=True), axis=1
    )
    median_ci_low, median_ci_high = np.quantile(bootstrap_medians, [0.025, 0.975])
    score_summary = pd.DataFrame([{
        "contrast": "paired_tumour_vs_normal", "n_pairs": len(wide_score),
        "median_delta": (wide_score["T"] - wide_score["N"]).median(),
        "median_delta_bootstrap_ci_low": median_ci_low,
        "median_delta_bootstrap_ci_high": median_ci_high,
        "mean_delta": (wide_score["T"] - wide_score["N"]).mean(),
        "n_higher_tumour": int((wide_score["T"] > wide_score["N"]).sum()),
        "wilcoxon_statistic": score_stat, "wilcoxon_p": score_p,
    }])
    score_summary.to_csv(OUT / "schwann_pseudobulk_paired_score_test.csv", index=False)

    gene_rows = []
    paired_long = []
    for gene in CORE:
        wide = pb.pivot(index="Donor", columns="Sample_type", values=f"logCPM_{gene}").dropna()
        delta = wide["T"] - wide["N"]
        stat, p = wilcoxon(wide["T"], wide["N"])
        gene_rows.append({
            "gene": gene, "n_pairs": len(wide), "median_delta_logCPM": delta.median(),
            "mean_delta_logCPM": delta.mean(), "n_higher_tumour": int((delta > 0).sum()),
            "wilcoxon_statistic": stat, "p": p,
        })
        paired_long.extend({"Donor": donor, "gene": gene, "delta_logCPM": value} for donor, value in delta.items())
    gene_tests = pd.DataFrame(gene_rows)
    gene_tests["fdr"] = multipletests(gene_tests["p"], method="fdr_bh")[1]
    gene_tests = gene_tests.sort_values("median_delta_logCPM", ascending=False)
    gene_tests.to_csv(OUT / "schwann_pseudobulk_gene_tests.csv", index=False)
    paired_long = pd.DataFrame(paired_long)

    tumour = pb[pb["Sample_type"].eq("T")]
    dmmr = tumour.loc[tumour["MMR status"].eq("dMMR"), "reactive_core_score"]
    pmmr = tumour.loc[tumour["MMR status"].eq("pMMR"), "reactive_core_score"]
    mmr_stat, mmr_p = mannwhitneyu(dmmr, pmmr, alternative="two-sided")
    pd.DataFrame([{
        "contrast": "dMMR_vs_pMMR_tumours_exploratory", "n_dMMR": len(dmmr), "n_pMMR": len(pmmr),
        "median_dMMR": dmmr.median(), "median_pMMR": pmmr.median(),
        "mannwhitney_statistic": mmr_stat, "p": mmr_p,
    }]).to_csv(OUT / "schwann_pseudobulk_mmr_test.csv", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1, 1.8]})
    plot_pair = wide_score.reset_index().melt(id_vars="Donor", value_vars=["N", "T"], var_name="Tissue", value_name="Score")
    for _, x in plot_pair.groupby("Donor"):
        axes[0].plot([0, 1], x.set_index("Tissue").loc[["N", "T"], "Score"], color="#8c8c8c", lw=1, alpha=.8)
    sns.stripplot(data=plot_pair, x="Tissue", y="Score", order=["N", "T"], hue="Tissue", palette=["#4C78A8", "#E45756"], size=7, ax=axes[0], legend=False)
    axes[0].set_xticks([0, 1], ["Normal", "Tumour"])
    axes[0].set_title(f"Paired donors (n={len(wide_score)})\nWilcoxon P={score_p:.3g}")
    axes[0].set_ylabel("13-gene scaffold pseudobulk score")
    axes[0].set_xlabel("")

    order = gene_tests["gene"].tolist()
    sns.stripplot(data=paired_long, x="delta_logCPM", y="gene", order=order, color="#777777", alpha=.55, size=4, ax=axes[1])
    med = paired_long.groupby("gene")["delta_logCPM"].median().reindex(order)
    axes[1].scatter(med, np.arange(len(order)), color="#D62728", s=42, zorder=5, label="Median")
    axes[1].axvline(0, color="black", lw=.8, ls="--")
    axes[1].set_xlabel("Tumour − normal logCPM (each point = donor)")
    axes[1].set_ylabel("")
    axes[1].set_title("Component-gene concordance")
    axes[1].legend(frameon=False)
    fig.suptitle("GSE303070 CosMx: tumour-associated reactive Schwann program", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "GSE303070_schwann_pseudobulk_validation.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "GSE303070_schwann_pseudobulk_validation.pdf", bbox_inches="tight")
    plt.close(fig)

    print(score_summary.to_string(index=False))
    print("\nPer-gene paired tests\n", gene_tests.to_string(index=False))
    print(f"\nExploratory MMR: dMMR n={len(dmmr)}, pMMR n={len(pmmr)}, P={mmr_p:.4g}")


if __name__ == "__main__":
    main()
