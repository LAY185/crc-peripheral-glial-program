#!/usr/bin/env python3
"""Assemble four submission-grade manuscript figures from locked result tables."""

from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
OUT = ROOT / "results/manuscript_figures"
SOURCE = ROOT / "results/source_data"

BLUE = "#3B6EA8"
RED = "#D35D4F"
GREY = "#777777"
LIGHT_GREY = "#D9D9D9"
GOLD = "#D9A441"
TEAL = "#3C8D8D"

CORE = [
    "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
    "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
]
MAPPING_PATH = ROOT / "config/mouse_human_scaffold_mapping.csv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})
sns.set_style("ticks")


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(
        OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def score_table(dataset: str, threshold: int = 3) -> pd.DataFrame:
    d = pd.read_csv(TABLES / f"{dataset}_marker_summary.csv")
    d = d[
        d["Cell_subtype"].eq("Enteric glial cells")
        & d["Class"].isin(["Normal", "Tumor"])
        & (d["n_cells"] >= threshold)
    ].copy()
    cols = [f"{gene}__mean" for gene in CORE if f"{gene}__mean" in d]
    d["score"] = d[cols].mean(axis=1)
    return d[["Patient", "Class", "n_cells", "score"]]


def figure1() -> None:
    mapping = pd.read_csv(MAPPING_PATH)
    required = {"mouse_symbol", "human_symbol", "ensembl_human_gene_id", "ensembl_orthology_type", "selection_note"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"Orthologue mapping lacks required columns: {sorted(required - set(mapping.columns))}")
    if mapping["mouse_symbol"].duplicated().any() or mapping["human_symbol"].duplicated().any():
        raise ValueError("Orthologue mapping contains duplicated mouse or human symbols")
    if set(mapping["human_symbol"]) != set(CORE):
        raise ValueError("Orthologue mapping and 13-gene scaffold do not match")
    if mapping.loc[mapping["human_symbol"].eq("MT2A"), "ensembl_human_gene_id"].item() != "ENSG00000125148":
        raise ValueError("MT2A Ensembl identifier mismatch")
    orthologue_map = mapping.set_index("mouse_symbol")["human_symbol"].to_dict()
    d = pd.read_csv(TABLES / "GSE231802_limma/reactive_egc_consensus.csv")
    d = d[d["n_up_fdr05"] >= 2].copy().sort_values("mean_logFC", ascending=False)
    d["mouse_symbol"] = d["gene"]
    d["human_orthologue"] = d["gene"].map(orthologue_map)
    d["cross_platform_core"] = d["human_orthologue"].notna()
    SOURCE.mkdir(parents=True, exist_ok=True)
    d.to_csv(SOURCE / "Figure1_reactive_program.csv", index=False)
    mapping.to_csv(SOURCE / "Figure1_mouse_human_orthologue_mapping.csv", index=False)

    fig = plt.figure(figsize=(7.20, 6.05), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 3.4], width_ratios=[2.2, 1.0])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(0, 2)
    ax_a.axis("off")
    items = [
        (0.2, "Mouse enteric glia", LIGHT_GREY),
        (2.45, "Tumour CM\nHealthy CM\nControl", "#E8EEF5"),
        (5.0, "6 h   12 h   24 h", "#F5EBDD"),
        (7.35, "limma contrasts", "#F4DEDA"),
    ]
    for x, text, color in items:
        box = FancyBboxPatch((x, 0.45), 1.75, 1.0, boxstyle="round,pad=0.04,rounding_size=0.05", fc=color, ec="#555555", lw=0.7)
        ax_a.add_patch(box)
        ax_a.text(x + 0.875, 0.95, text, ha="center", va="center", fontsize=7)
    for x in [2.05, 4.3, 6.85]:
        ax_a.annotate("", xy=(x + 0.28, 0.95), xytext=(x, 0.95), arrowprops={"arrowstyle": "->", "lw": 0.8, "color": GREY})
    ax_a.text(9.48, 0.95, "33-gene\nprogram", ha="center", va="center", fontsize=7, fontweight="bold", color=RED)
    panel_label(ax_a, "a")

    heat = d.set_index("mouse_symbol")[["logFC_6h", "logFC_12h", "logFC_24h"]]
    heat.columns = ["6 h", "12 h", "24 h"]
    norm = TwoSlopeNorm(vmin=min(-1, float(heat.min().min())), vcenter=0, vmax=float(heat.max().max()))
    im = ax_b.imshow(heat.to_numpy(), aspect="auto", cmap="RdBu_r", norm=norm)
    ax_b.set_xticks(range(3), heat.columns)
    ax_b.set_yticks(range(len(heat)), heat.index, fontsize=5.4)
    ax_b.set_ylabel("Mouse gene symbol")
    ax_b.tick_params(length=0)
    ax_b.set_title("Time-consistent tumour-conditioned response", loc="left", pad=5)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.035, pad=0.02)
    cbar.set_label("log2 fold change", fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    panel_label(ax_b, "b")

    ax_c.axis("off")
    ax_c.set_title("Fixed cross-platform scaffold", loc="left", pad=5)
    ax_c.text(0.03, 0.96, "33 perturbation-derived genes", transform=ax_c.transAxes, fontsize=7, va="top")
    ax_c.annotate("", xy=(0.50, 0.84), xytext=(0.50, 0.91), xycoords="axes fraction", arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.8})
    ax_c.text(0.03, 0.82, "Restrict to genes measured by\nboth CosMx panel versions", transform=ax_c.transAxes, fontsize=6.5, va="top")
    ax_c.annotate("", xy=(0.50, 0.67), xytext=(0.50, 0.74), xycoords="axes fraction", arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.8})
    core_text = "\n".join([f"{i + 1:>2}. {g}" for i, g in enumerate(CORE)])
    ax_c.text(0.08, 0.64, core_text, transform=ax_c.transAxes, family="monospace", fontsize=6.2, va="top", color="#222222")
    ax_c.text(0.03, 0.05, "Selection precedes all human comparisons", transform=ax_c.transAxes, fontsize=6.2, color=RED)
    panel_label(ax_c, "c")
    save_figure(fig, "Figure1_reactive_program_derivation")


def draw_patient_scores(ax, d: pd.DataFrame, dataset: str) -> None:
    mapping = {"Normal": 0, "Tumor": 1}
    for patient, x in d.groupby("Patient"):
        if set(x["Class"]) == {"Normal", "Tumor"}:
            y = x.set_index("Class").loc[["Normal", "Tumor"], "score"]
            ax.plot([0, 1], y, color=LIGHT_GREY, lw=0.7, zorder=1)
    for cls, color in [("Normal", BLUE), ("Tumor", RED)]:
        x = d[d["Class"].eq(cls)]
        jitter = np.linspace(-0.06, 0.06, len(x)) if len(x) > 1 else np.array([0.0])
        sizes = 12 + 16 * np.sqrt(x["n_cells"] / x["n_cells"].max())
        ax.scatter(mapping[cls] + jitter, x["score"], s=sizes, c=color, ec="white", lw=0.4, alpha=0.9, zorder=2)
        ax.hlines(x["score"].median(), mapping[cls] - 0.18, mapping[cls] + 0.18, color="#222222", lw=1.1, zorder=3)
    ax.set_xticks([0, 1], ["Normal", "Tumour"])
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel("13-gene scaffold score")
    ax.set_title(dataset, loc="left")
    normal_median = d.loc[d["Class"].eq("Normal"), "score"].median()
    tumour_median = d.loc[d["Class"].eq("Tumor"), "score"].median()
    ax.text(
        0.03, 0.96,
        f"Normal median = {normal_median:.3f}\nTumour median = {tumour_median:.3f}\nMedian difference = {tumour_median - normal_median:.3f}",
        transform=ax.transAxes, va="top", fontsize=6,
    )


def figure2() -> None:
    atlas = TABLES / "crc_atlas_external_validation"
    deltas = pd.read_csv(atlas / "paired_patient_deltas.csv")
    tests = pd.read_csv(atlas / "paired_patient_tests.csv")
    studies = pd.read_csv(atlas / "study_level_summary.csv")
    loo = pd.read_csv(atlas / "leave_one_study_out.csv")
    primary = deltas[
        deltas["min_cells_per_sample"].eq(5)
        & deltas["score"].eq("three_gene")
        & deltas["model"].eq("depth_adjusted")
    ].copy()
    main_test = tests[
        tests["min_cells_per_sample"].eq(5)
        & tests["score"].eq("three_gene")
        & tests["model"].eq("depth_adjusted")
    ].iloc[0]
    study = studies[studies["score"].eq("three_gene")].copy().sort_values("median_delta")
    sensitivity = tests[
        tests["score"].eq("three_gene") & tests["model"].eq("depth_adjusted")
    ].copy().sort_values("min_cells_per_sample")
    loo3 = loo[loo["score"].eq("three_gene")].copy().sort_values("median_delta")
    SOURCE.mkdir(parents=True, exist_ok=True)
    primary.to_csv(SOURCE / "Figure2_CRC_atlas_paired_patients.csv", index=False)
    study.to_csv(SOURCE / "Figure2_CRC_atlas_studies.csv", index=False)
    sensitivity.to_csv(SOURCE / "Figure2_CRC_atlas_threshold_sensitivity.csv", index=False)
    loo3.to_csv(SOURCE / "Figure2_CRC_atlas_leave_one_study_out.csv", index=False)

    fig = plt.figure(figsize=(7.20, 5.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.35], height_ratios=[1.1, 0.9])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    for _, row in primary.iterrows():
        ax_a.plot([0, 1], [row["normal"], row["tumour"]], color=LIGHT_GREY, lw=0.55, zorder=1)
    ax_a.scatter(np.zeros(len(primary)), primary["normal"], c=BLUE, s=14, ec="white", lw=0.3, zorder=2)
    ax_a.scatter(np.ones(len(primary)), primary["tumour"], c=RED, s=14, ec="white", lw=0.3, zorder=2)
    ax_a.set_xticks([0, 1], ["Adjacent normal", "Primary tumour"])
    ax_a.set_xlim(-0.35, 1.35)
    ax_a.set_ylabel("Depth-adjusted three-gene score")
    ax_a.set_title("Paired patients across CRC Atlas", loc="left")
    ax_a.text(
        0.03, 0.98,
        "38 of 44 patients increase\n9 independent studies\nMedian Δ = 0.124\nHierarchical bootstrap 95% CI 0.035 to 0.198\nWilcoxon P = 4.60e-7",
        transform=ax_a.transAxes, va="top", fontsize=6,
    )
    panel_label(ax_a, "a")

    y = np.arange(len(study))
    ax_b.axvline(0, color=GREY, lw=0.7, ls=":")
    ax_b.scatter(study["median_delta"], y, s=18 + 2.0 * study["n_paired_patients"], color=TEAL, ec="white", lw=0.4)
    ax_b.set_yticks(y, [x.replace("_", " ") for x in study["study_id"]])
    ax_b.set_xlabel("Median paired tumour-normal Δ")
    ax_b.set_title("All study-level medians are positive", loc="left")
    panel_label(ax_b, "b")

    ax = ax_c
    x = sensitivity["min_cells_per_sample"].to_numpy(float)
    y = sensitivity["median_delta"].to_numpy(float)
    lo = sensitivity["hierarchical_bootstrap_ci_low"].to_numpy(float)
    hi = sensitivity["hierarchical_bootstrap_ci_high"].to_numpy(float)
    ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", ms=4, color=RED, ecolor=RED, capsize=2, lw=0.9)
    ax.axhline(0, color=GREY, lw=0.7, ls=":")
    ax.set_xticks(x)
    ax.set_xlabel("Minimum Schwann cells per tissue sample")
    ax.set_ylabel("Median paired Δ (95% hierarchical bootstrap CI)")
    ax.set_title("Cell-count threshold sensitivity", loc="left")
    for _, row in sensitivity.iterrows():
        ax.text(row["min_cells_per_sample"], row["hierarchical_bootstrap_ci_high"] + 0.012, f"{int(row['n_positive'])}/{int(row['n_paired_patients'])}", ha="center", va="bottom", fontsize=5.3, color=GREY)
    panel_label(ax, "c")
    ax.text(0.98, 0.03, f"Leave-one-study-out median Δ {loo3['median_delta'].min():.3f} to {loo3['median_delta'].max():.3f}\nP ≤ {loo3['p'].max():.2g}", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.3, color=GREY)
    save_figure(fig, "Figure2_crc_atlas_validation")


def figure3() -> None:
    pb = pd.read_csv(TABLES / "GSE303070_cosmx/schwann_pseudobulk_scores.csv")
    tests = pd.read_csv(TABLES / "GSE303070_cosmx/schwann_pseudobulk_gene_tests.csv")
    wide_all = pb.pivot(index="Donor", columns="Sample_type", values="reactive_core_score")
    wide = wide_all.loc[wide_all[["N", "T"]].notna().all(axis=1)]
    delta_columns = {}
    for gene in CORE:
        gene_wide = pb.pivot(index="Donor", columns="Sample_type", values=f"logCPM_{gene}")
        gene_wide = gene_wide.loc[gene_wide[["N", "T"]].notna().all(axis=1)]
        delta_columns[gene] = gene_wide["T"] - gene_wide["N"]
    delta = pd.DataFrame(delta_columns)
    order = tests.sort_values("median_delta_logCPM", ascending=False)["gene"].tolist()
    delta = delta[order]
    SOURCE.mkdir(parents=True, exist_ok=True)
    wide.reset_index().to_csv(SOURCE / "Figure3_paired_scores.csv", index=False)
    delta.reset_index().to_csv(SOURCE / "Figure3_paired_gene_deltas.csv", index=False)
    tests.to_csv(SOURCE / "Figure3_gene_tests.csv", index=False)
    nomination = tests.copy()
    nomination["positive_q_rank"] = np.nan
    positive = nomination[nomination["median_delta_logCPM"].gt(0)].sort_values(
        ["fdr", "n_higher_tumour", "median_delta_logCPM"],
        ascending=[True, False, False],
    )
    nomination.loc[positive.index, "positive_q_rank"] = np.arange(1, len(positive) + 1)
    minimum_positive_q = positive["fdr"].min()
    nomination["selected_for_atlas"] = (
        nomination["median_delta_logCPM"].gt(0)
        & np.isclose(nomination["fdr"], minimum_positive_q)
    )
    selected = nomination[nomination["selected_for_atlas"]].sort_values("positive_q_rank")["gene"].tolist()
    if selected != ["TIMP1", "IFITM1", "CHI3L1"]:
        raise ValueError(f"Unexpected CosMx nomination: {selected}")
    nomination = nomination.rename(columns={"p": "wilcoxon_p", "fdr": "bh_q"})
    nomination = nomination[
        ["gene", "n_pairs", "median_delta_logCPM", "n_higher_tumour", "wilcoxon_p", "bh_q", "positive_q_rank", "selected_for_atlas"]
    ].sort_values(["selected_for_atlas", "positive_q_rank", "median_delta_logCPM"], ascending=[False, True, False])
    nomination.to_csv(SOURCE / "Figure3_cosmx_component_nomination.csv", index=False)

    fig = plt.figure(figsize=(7.20, 4.45), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.9, 2.1], height_ratios=[1, 1])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    for _, row in wide.iterrows():
        ax_a.plot([0, 1], [row["N"], row["T"]], color=LIGHT_GREY, lw=0.8, zorder=1)
    ax_a.scatter(np.zeros(len(wide)), wide["N"], c=BLUE, s=25, ec="white", lw=0.4, zorder=2)
    ax_a.scatter(np.ones(len(wide)), wide["T"], c=RED, s=25, ec="white", lw=0.4, zorder=2)
    ax_a.set_xticks([0, 1], ["Normal", "Tumour"])
    ax_a.set_xlim(-0.35, 1.35)
    ax_a.set_ylabel("Mean 13-gene scaffold log1p CPM")
    ax_a.set_title("Paired Schwann pseudobulk", loc="left")
    ax_a.text(0.03, 0.98, "7 of 8 donors increase\nMedian Δ = 0.168\nBootstrap 95% interval 0.100 to 0.506\nWilcoxon P = 0.0234", transform=ax_a.transAxes, va="top", fontsize=6)
    panel_label(ax_a, "a")

    vmax = float(np.abs(delta.to_numpy()).max())
    im = ax_b.imshow(delta.to_numpy(), aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax_b.set_xticks(range(len(order)), order, rotation=55, ha="right", rotation_mode="anchor")
    ax_b.set_yticks(range(len(delta)), [x.replace("Donor", "D") for x in delta.index])
    ax_b.set_title("Tumour-normal change within each donor", loc="left")
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.025, pad=0.02)
    cbar.set_label("Δ log1p CPM", fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    panel_label(ax_b, "b")

    plot = tests.set_index("gene").loc[order].iloc[::-1]
    y = np.arange(len(plot))
    colors = [RED if v > 0 else BLUE for v in plot["median_delta_logCPM"]]
    ax_c.axvline(0, color=GREY, lw=0.7, ls=":")
    ax_c.hlines(y, 0, plot["median_delta_logCPM"], color=LIGHT_GREY, lw=1)
    ax_c.scatter(plot["median_delta_logCPM"], y, c=colors, s=20, zorder=2)
    ax_c.set_yticks(y, plot.index)
    ax_c.set_xlabel("Median paired Δ log1p CPM")
    ax_c.set_title("All 13 component genes", loc="left")
    for yi, (_, row) in enumerate(plot.iterrows()):
        ax_c.text(ax_c.get_xlim()[1], yi, f"q={row['fdr']:.3f}", va="center", ha="right", fontsize=5.2, color=GREY)
    panel_label(ax_c, "c")
    save_figure(fig, "Figure3_cosmx_paired_validation")


def figure4() -> None:
    a = score_table("GSE132465")
    b = score_table("GSE144735")
    sens = pd.read_csv(TABLES / "cross_platform_meta/scrna_cell_threshold_sensitivity.csv")
    spatial = pd.read_csv(TABLES / "commsbio_spatial_gradient/specimen_spatial_coefficients.csv")
    spatial = spatial[
        spatial["score"].eq("core13")
        & spatial["proximity"].eq("epi_k20")
        & spatial["Sample_type"].eq("T")
    ].copy()
    SOURCE.mkdir(parents=True, exist_ok=True)
    a.assign(dataset="GSE132465").to_csv(SOURCE / "Figure4_GSE132465_scores.csv", index=False)
    b.assign(dataset="GSE144735").to_csv(SOURCE / "Figure4_GSE144735_scores.csv", index=False)
    sens.to_csv(SOURCE / "Figure4_scrna_threshold_sensitivity.csv", index=False)
    spatial.to_csv(SOURCE / "Figure4_spatial_boundary.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(7.20, 5.15), constrained_layout=True)
    draw_patient_scores(axes[0, 0], a, "GSE132465")
    draw_patient_scores(axes[0, 1], b, "GSE144735")
    axes[0, 1].set_ylabel("")
    panel_label(axes[0, 0], "a")
    panel_label(axes[0, 1], "b")

    ax = axes[1, 0]
    offsets = {"GSE132465": -0.08, "GSE144735": 0.08}
    colors = {"GSE132465": BLUE, "GSE144735": RED}
    for dataset, xdf in sens.groupby("dataset"):
        valid = xdf["hedges_g"].notna()
        xv = xdf.loc[valid, "minimum_glial_cells"].to_numpy(float) + offsets[dataset]
        yv = xdf.loc[valid, "hedges_g"].to_numpy(float)
        ax.scatter(xv, yv, s=18, color=colors[dataset], label=dataset)
    ax.axhline(0, color=GREY, lw=0.7, ls=":")
    ax.set_xticks([1, 3, 5, 10])
    ax.set_xlabel("Minimum enteric-glial cells per summary")
    ax.set_ylabel("Descriptive Hedges g")
    ax.set_title("Positive standardized differences across thresholds", loc="left")
    ax.legend(loc="upper right", frameon=False)
    panel_label(ax, "c")

    ax = axes[1, 1]
    spatial_summary = []
    rng = np.random.default_rng(20260906)
    for i, (model, color) in enumerate([("unadjusted", GOLD), ("adjusted", TEAL)]):
        vals = spatial.loc[spatial["model"].eq(model), "standardized_beta"]
        boot = np.median(rng.choice(vals.to_numpy(), size=(200_000, len(vals)), replace=True), axis=1)
        ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
        spatial_summary.append({
            "model": model, "n_specimens": len(vals), "median_standardized_beta": vals.median(),
            "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
        })
        jitter = np.linspace(-0.055, 0.055, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=15, color=color, alpha=0.7, ec="white", lw=0.3)
        ax.scatter(i, vals.median(), s=42, marker="D", color=color, ec="white", lw=0.5, zorder=3)
    ax.axhline(0, color=GREY, lw=0.7, ls=":")
    ax.set_xticks([0, 1], ["Unadjusted", "Adjusted"])
    ax.set_ylabel("Specimen-level standardized β")
    ax.set_title("Epithelial proximity loses support after adjustment", loc="left")
    pd.DataFrame(spatial_summary).to_csv(SOURCE / "Figure4_spatial_boundary_summary.csv", index=False)
    ax.text(0.03, 0.97, "Unadjusted median β = 0.116\nAdjusted median β = 0.00180\nAdjusted bootstrap 95% CI -0.048 to 0.071\nFOV permutation P = 0.876", transform=ax.transAxes, va="top", fontsize=5.3)
    panel_label(ax, "d")
    save_figure(fig, "Figure4_supportive_cohorts_and_boundary")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    figure1()
    figure2()
    figure3()
    figure4()
    shutil.copy2(ROOT / "config/figure_legends.md", OUT / "figure_legends.md")
    print(f"Wrote manuscript figures to {OUT}")


if __name__ == "__main__":
    main()
