#!/usr/bin/env python3
"""Independent patient-level validation in the 4.27M-cell CRC Atlas."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/processed/CRC_atlas_schwann_core13_expression.csv.gz"
OUT = ROOT / "results/tables/crc_atlas_external_validation"
FIG = ROOT / "results/figures"
EXCLUDED_OVERLAP = {"Lee_2020_Nat_Genet"}
GENESETS = {
    "core13": [
        "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
        "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
    ],
    "noncytokine8": ["MT2A", "SOD2", "IER3", "LCN2", "SLPI", "CHI3L1", "TIMP1", "IFITM1"],
    "three_gene": ["TIMP1", "IFITM1", "CHI3L1"],
}
THRESHOLDS = [3, 5, 10, 20]
RNG = np.random.default_rng(20260904)


def depth_residualize(frame: pd.DataFrame, score_col: str) -> pd.Series:
    result = pd.Series(index=frame.index, dtype=float)
    for _, idx in frame.groupby("study_id", observed=True).groups.items():
        part = frame.loc[idx]
        x = np.column_stack([
            np.ones(len(part)),
            np.log1p(part["n_counts"].to_numpy(float)),
            np.log1p(part["n_genes"].to_numpy(float)),
        ])
        y = part[score_col].to_numpy(float)
        coef = np.linalg.lstsq(x, y, rcond=None)[0]
        result.loc[idx] = y - x @ coef + np.mean(y)
    return result


def paired_deltas(samples: pd.DataFrame, threshold: int, value_col: str) -> pd.DataFrame:
    eligible = samples[samples["n_cells"] >= threshold]
    patient = (
        eligible.groupby(["study_id", "patient_id", "sample_type"], observed=True)[value_col]
        .mean().reset_index()
    )
    patient = patient[patient["sample_type"].isin(["primary tumor", "adjacent normal"])]
    wide = patient.pivot(index=["study_id", "patient_id"], columns="sample_type", values=value_col)
    required = ["primary tumor", "adjacent normal"]
    if not set(required).issubset(wide.columns):
        return pd.DataFrame(columns=["study_id", "patient_id", "normal", "tumour", "delta"])
    wide = wide.dropna(subset=required).reset_index()
    return wide.rename(columns={"adjacent normal": "normal", "primary tumor": "tumour"}).assign(
        delta=lambda x: x["tumour"] - x["normal"]
    )


def hierarchical_bootstrap_ci(deltas: pd.DataFrame, n_boot: int = 5000) -> tuple[float, float]:
    studies = deltas["study_id"].unique()
    boot = np.empty(n_boot)
    grouped = {study: x["delta"].to_numpy() for study, x in deltas.groupby("study_id")}
    for i in range(n_boot):
        picked = RNG.choice(studies, size=len(studies), replace=True)
        values = []
        for study in picked:
            source = grouped[study]
            values.extend(RNG.choice(source, size=len(source), replace=True))
        boot[i] = np.median(values)
    return tuple(np.quantile(boot, [0.025, 0.975]))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(INPUT)
    cells = cells[~cells["study_id"].isin(EXCLUDED_OVERLAP)].copy()

    score_names = []
    for name, genes in GENESETS.items():
        raw = f"{name}_raw"
        adjusted = f"{name}_depth_adjusted"
        cells[raw] = cells[genes].mean(axis=1)
        cells[adjusted] = depth_residualize(cells, raw)
        score_names.extend([raw, adjusted])

    samples = (
        cells.groupby(["study_id", "patient_id", "sample_id", "sample_type"], observed=True)
        .agg(n_cells=("obs_index", "size"), **{x: (x, "mean") for x in score_names})
        .reset_index()
    )
    samples.to_csv(OUT / "schwann_sample_scores.csv", index=False)

    test_rows, delta_frames = [], []
    for threshold in THRESHOLDS:
        for score in GENESETS:
            for model in ["raw", "depth_adjusted"]:
                value_col = f"{score}_{model}"
                delta = paired_deltas(samples, threshold, value_col)
                if delta.empty:
                    continue
                stat, p = wilcoxon(delta["delta"], alternative="two-sided")
                lo, hi = hierarchical_bootstrap_ci(delta)
                test_rows.append({
                    "min_cells_per_sample": threshold,
                    "score": score,
                    "model": model,
                    "n_paired_patients": len(delta),
                    "n_studies": delta["study_id"].nunique(),
                    "n_positive": int((delta["delta"] > 0).sum()),
                    "median_delta": delta["delta"].median(),
                    "hierarchical_bootstrap_ci_low": lo,
                    "hierarchical_bootstrap_ci_high": hi,
                    "wilcoxon_statistic": stat,
                    "p": p,
                })
                tagged = delta.assign(min_cells_per_sample=threshold, score=score, model=model)
                delta_frames.append(tagged)
    tests = pd.DataFrame(test_rows)
    tests["q_within_threshold_model"] = np.nan
    for _, idx in tests.groupby(["min_cells_per_sample", "model"]).groups.items():
        tests.loc[idx, "q_within_threshold_model"] = multipletests(tests.loc[idx, "p"], method="fdr_bh")[1]
    tests.to_csv(OUT / "paired_patient_tests.csv", index=False)
    deltas_all = pd.concat(delta_frames, ignore_index=True)
    deltas_all.to_csv(OUT / "paired_patient_deltas.csv", index=False)

    primary = deltas_all[
        (deltas_all["min_cells_per_sample"] == 5)
        & (deltas_all["model"] == "depth_adjusted")
    ].copy()
    study_summary = (
        primary.groupby(["score", "study_id"], observed=True)["delta"]
        .agg(n_paired_patients="size", median_delta="median", n_positive=lambda x: int((x > 0).sum()))
        .reset_index()
    )
    study_summary.to_csv(OUT / "study_level_summary.csv", index=False)

    loo_rows = []
    for score, score_data in primary.groupby("score", observed=True):
        for omitted in sorted(score_data["study_id"].unique()):
            remain = score_data[score_data["study_id"] != omitted]
            stat, p = wilcoxon(remain["delta"], alternative="two-sided")
            loo_rows.append({
                "score": score,
                "omitted_study": omitted,
                "n_paired_patients": len(remain),
                "n_studies": remain["study_id"].nunique(),
                "n_positive": int((remain["delta"] > 0).sum()),
                "median_delta": remain["delta"].median(),
                "p": p,
            })
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(OUT / "leave_one_study_out.csv", index=False)

    plot = primary.copy()
    order = (
        plot[plot["score"] == "three_gene"].groupby("study_id")["delta"].median()
        .sort_values().index.tolist()
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))

    three = plot[plot["score"] == "three_gene"]
    sns.stripplot(data=three, y="study_id", x="delta", order=order, ax=axes[0, 0],
                  color="#4C78A8", size=5, jitter=.16)
    axes[0, 0].axvline(0, color="black", lw=.8, ls="--")
    axes[0, 0].set(title="a  Independent paired patients", xlabel="Adjusted 3-gene tumour−normal score", ylabel="Study")

    st = study_summary[study_summary["score"] == "three_gene"].set_index("study_id").loc[order].reset_index()
    colors = np.where(st["median_delta"] > 0, "#E45756", "#9D9D9D")
    axes[0, 1].barh(st["study_id"], st["median_delta"], color=colors)
    axes[0, 1].axvline(0, color="black", lw=.8)
    for y, row in st.iterrows():
        axes[0, 1].text(row["median_delta"], y, f"  {int(row.n_positive)}/{int(row.n_paired_patients)}",
                        va="center", ha="left" if row["median_delta"] >= 0 else "right", fontsize=8)
    axes[0, 1].set(title="b  Study-level consistency", xlabel="Median adjusted paired difference", ylabel="")

    t5 = tests[tests["min_cells_per_sample"] == 5].copy()
    t5["label"] = t5["score"].map({"core13": "13-gene core", "noncytokine8": "Non-cytokine 8", "three_gene": "TIMP1/IFITM1/CHI3L1"})
    for i, model in enumerate(["raw", "depth_adjusted"]):
        x = t5[t5["model"] == model].set_index("score").loc[list(GENESETS)].reset_index()
        ypos = np.arange(3) + (i - .5) * .16
        axes[1, 0].errorbar(x["median_delta"], ypos,
                            xerr=[x["median_delta"] - x["hierarchical_bootstrap_ci_low"],
                                  x["hierarchical_bootstrap_ci_high"] - x["median_delta"]],
                            fmt="o", capsize=3, label=model.replace("_", " "))
    axes[1, 0].axvline(0, color="black", lw=.8, ls="--")
    axes[1, 0].set_yticks(np.arange(3), t5[t5["model"] == "raw"].set_index("score").loc[list(GENESETS), "label"])
    axes[1, 0].set(title="c  Prespecified score sensitivity", xlabel="Median paired difference (95% hierarchical bootstrap CI)", ylabel="")
    axes[1, 0].legend(frameon=False, fontsize=8)

    sens = tests[(tests["score"] == "three_gene") & (tests["model"] == "depth_adjusted")]
    axes[1, 1].plot(sens["min_cells_per_sample"], sens["median_delta"], marker="o", color="#E45756")
    axes[1, 1].fill_between(sens["min_cells_per_sample"], sens["hierarchical_bootstrap_ci_low"],
                            sens["hierarchical_bootstrap_ci_high"], color="#E45756", alpha=.18)
    for _, row in sens.iterrows():
        offsets = {3: (-18, 18), 5: (18, 5), 10: (0, 8), 20: (0, 8)}
        axes[1, 1].annotate(f"n={int(row.n_paired_patients)}\nP={row.p:.2g}",
                            (row.min_cells_per_sample, row.median_delta),
                            xytext=offsets[int(row.min_cells_per_sample)],
                            textcoords="offset points", ha="center", fontsize=8)
    axes[1, 1].axhline(0, color="black", lw=.8, ls="--")
    axes[1, 1].set(title="d  Cell-count threshold sensitivity", xlabel="Minimum Schwann cells per sample", ylabel="Median adjusted paired difference")

    sns.despine(fig)
    fig.suptitle("CRC Atlas external validation (Lee 2020 overlap excluded)", y=.995, fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG / "CommsBio_gate_crc_atlas_validation.png", dpi=300)
    fig.savefig(FIG / "CommsBio_gate_crc_atlas_validation.pdf")
    plt.close(fig)

    print(tests[tests["min_cells_per_sample"] == 5].to_string(index=False))
    print("\nThree-gene leave-one-study-out P range:",
          loo.loc[loo["score"] == "three_gene", "p"].min(),
          loo.loc[loo["score"] == "three_gene", "p"].max())


if __name__ == "__main__":
    main()
