#!/usr/bin/env python3
"""Communications Biology gate 2: tumour-proximity spatial gradients.

For each author-annotated Schwann cell, derive continuous spatial exposure
metrics inside its CosMx FOV. Estimate one within-FOV adjusted coefficient per
biological specimen, then test those specimen coefficients. Cell-level P values
are deliberately not used.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data/raw/GSE303070/All_cells_metadata.csv.gz"
SCHWANN = ROOT / "data/processed/GSE303070_schwann_identity_targets.csv.gz"
OUT = ROOT / "results/tables/commsbio_spatial_gradient"
FIG = ROOT / "results/figures"
K = 20
N_PERM = 2000
CORE = ["MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI", "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1"]
NO_CYTOKINES = [g for g in CORE if g not in {"CCL2", "CXCL10", "IL6", "CSF3", "CCL11"}]
THREE_GENE = ["TIMP1", "IFITM1", "CHI3L1"]
MYELOID = ["PTPRC", "LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "CD68", "CD163", "APOE"]
FIBROBLAST = ["COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"]


def score_available(d: pd.DataFrame, genes: list[str]) -> pd.Series:
    cols = [f"expr_{g}" for g in genes if f"expr_{g}" in d and d[f"expr_{g}"].notna().any()]
    return d[cols].mean(axis=1, skipna=True)


def add_spatial_exposures(meta: pd.DataFrame, schwann_barcodes: set[str]) -> pd.DataFrame:
    records = []
    for (slide, fov), x in meta.groupby(["Slide", "fov"], sort=False):
        focal_mask = x["Barcode"].isin(schwann_barcodes)
        if not focal_mask.any() or len(x) <= K:
            continue
        x = x.copy()
        xc = x["CenterX_local_px"].fillna(x["CenterX_global_px"])
        yc = x["CenterY_local_px"].fillna(x["CenterY_global_px"])
        valid = xc.notna() & yc.notna()
        x = x.loc[valid].copy()
        if len(x) <= K:
            continue
        xy = np.column_stack([xc.loc[valid].to_numpy(float), yc.loc[valid].to_numpy(float)])
        focal_idx = np.flatnonzero(x["Barcode"].isin(schwann_barcodes).to_numpy())
        if not len(focal_idx):
            continue
        tree = cKDTree(xy)
        nn = tree.query(xy[focal_idx], k=min(K + 1, len(x)))[1]
        # Remove the focal cell itself, which is always the nearest point.
        nn = nn[:, 1:]
        labels = x["Manual_toplevel_pred"].astype(str).to_numpy()
        niches = x["niches"].fillna("NA").astype(str).to_numpy()
        neigh_labels = labels[nn]
        neigh_niches = niches[nn]

        epi_idx = np.flatnonzero(labels == "Epi")
        epi_mass_idx = np.flatnonzero(niches == "Epithelial mass")
        epi_dist = cKDTree(xy[epi_idx]).query(xy[focal_idx], k=1)[0] if len(epi_idx) else np.full(len(focal_idx), np.nan)
        mass_dist = cKDTree(xy[epi_mass_idx]).query(xy[focal_idx], k=1)[0] if len(epi_mass_idx) else np.full(len(focal_idx), np.nan)

        rec = pd.DataFrame({
            "Barcode": x.iloc[focal_idx]["Barcode"].to_numpy(),
            "Slide": slide,
            "fov": fov,
            "epi_k20_fraction": (neigh_labels == "Epi").mean(axis=1),
            "myeloid_k20_fraction": np.isin(neigh_labels, ["Macro", "Mono", "DC"]).mean(axis=1),
            "fibro_k20_fraction": (neigh_labels == "Fibro").mean(axis=1),
            "epithelial_mass_k20_fraction": (neigh_niches == "Epithelial mass").mean(axis=1),
            "nearest_epi_distance_px": epi_dist,
            "nearest_epithelial_mass_distance_px": mass_dist,
        })
        records.append(rec)
    return pd.concat(records, ignore_index=True)


def within_fov_demean(x: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    return x[cols] - x.groupby(["Slide", "fov"])[cols].transform("mean")


def standardize_frame(x: pd.DataFrame) -> pd.DataFrame:
    sd = x.std(axis=0, ddof=1).replace(0, np.nan)
    return (x - x.mean(axis=0)) / sd


def specimen_coefficient(x: pd.DataFrame, outcome: str, proximity: str, adjusted: bool) -> tuple[float, int, int]:
    nuisance = ["log_library", "expr_S100B", "myeloid_expr", "fibroblast_expr", "epithelial_expr", "log_panck"]
    cols = [outcome, proximity] + (nuisance if adjusted else [])
    y = x.dropna(subset=cols).copy()
    keep_fov = y.groupby(["Slide", "fov"])["Barcode"].transform("size").ge(5)
    y = y.loc[keep_fov]
    if len(y) < 30 or y[["Slide", "fov"]].drop_duplicates().shape[0] < 2:
        return np.nan, len(y), y[["Slide", "fov"]].drop_duplicates().shape[0]
    dm = within_fov_demean(y, cols)
    z = standardize_frame(dm).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 30 or z[proximity].std() == 0:
        return np.nan, len(z), y.loc[z.index, ["Slide", "fov"]].drop_duplicates().shape[0]
    predictors = [proximity] + (nuisance if adjusted else [])
    X = z[predictors].to_numpy(float)
    yy = z[outcome].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    return float(beta[0]), len(z), y.loc[z.index, ["Slide", "fov"]].drop_duplicates().shape[0]


def coefficient_table(d: pd.DataFrame) -> pd.DataFrame:
    score_map = {"core13": "score_core13", "noncytokine8": "score_noncytokine8", "three_gene": "score_three_gene"}
    proximity_map = {
        "epi_k20": "epi_k20_fraction",
        "nearest_epi": "negative_log_nearest_epi",
        "epithelial_mass_k20": "epithelial_mass_k20_fraction",
        "nearest_epithelial_mass": "negative_log_nearest_epithelial_mass",
    }
    rows = []
    for (sample, donor, tissue), x in d.groupby(["Sample", "Donor", "Sample_type"], sort=False):
        for score_name, outcome in score_map.items():
            for prox_name, proximity in proximity_map.items():
                for adjusted in [False, True]:
                    beta, n_cells, n_fovs = specimen_coefficient(x, outcome, proximity, adjusted)
                    rows.append({
                        "Sample": sample, "Donor": donor, "Sample_type": tissue,
                        "score": score_name, "proximity": prox_name,
                        "model": "adjusted" if adjusted else "unadjusted",
                        "standardized_beta": beta, "n_cells": n_cells, "n_fovs": n_fovs,
                    })
    return pd.DataFrame(rows)


def specimen_level_tests(coefs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, x in coefs.groupby(["score", "proximity", "model"]):
        score, proximity, model = keys
        for tissue in ["T", "N"]:
            vals = x.loc[x["Sample_type"].eq(tissue), "standardized_beta"].dropna()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                stat, p = wilcoxon(vals) if len(vals) else (np.nan, np.nan)
            rows.append({"score": score, "proximity": proximity, "model": model, "contrast": f"{tissue}_vs_zero", "n": len(vals), "median_beta": vals.median(), "wilcoxon_statistic": stat, "p": p})
        wide = x.pivot(index="Donor", columns="Sample_type", values="standardized_beta").dropna(subset=["N", "T"])
        delta = wide["T"] - wide["N"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stat, p = wilcoxon(delta) if len(delta) else (np.nan, np.nan)
        rows.append({"score": score, "proximity": proximity, "model": model, "contrast": "paired_T_minus_N", "n": len(delta), "median_beta": delta.median(), "wilcoxon_statistic": stat, "p": p})
    out = pd.DataFrame(rows)
    # Multiplicity is controlled over the four spatial metrics within each
    # score/model/contrast family. The three score definitions are sensitivity analyses.
    out["q_within_spatial_family"] = out.groupby(["score", "model", "contrast"])["p"].transform(
        lambda p: bh(p)
    )
    return out


def bh(p: pd.Series) -> pd.Series:
    arr = p.to_numpy(float)
    ok = np.isfinite(arr)
    ans = np.full(len(arr), np.nan)
    if ok.any():
        pv = arr[ok]
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(len(pv))
        tmp[order] = np.minimum(ranked, 1)
        ans[ok] = tmp
    return pd.Series(ans, index=p.index)


def permute_global_tumour(d: pd.DataFrame, outcome="score_core13", proximity="epi_k20_fraction"):
    tumour = d[d["Sample_type"].eq("T")].copy()
    observed_rows = []
    usable = {}
    for key, x in tumour.groupby(["Sample", "Donor"], sort=False):
        beta, n, nf = specimen_coefficient(x, outcome, proximity, adjusted=True)
        if np.isfinite(beta):
            observed_rows.append({"Sample": key[0], "Donor": key[1], "beta": beta, "n_cells": n, "n_fovs": nf})
            usable[key] = x.copy()
    obs = pd.DataFrame(observed_rows)
    observed_median = obs["beta"].median()
    rng = np.random.default_rng(20260904)
    null = np.empty(N_PERM)
    for b in range(N_PERM):
        betas = []
        for x in usable.values():
            xp = x.copy()
            xp[outcome] = xp.groupby(["Slide", "fov"])[outcome].transform(
                lambda v: rng.permutation(v.to_numpy())
            )
            beta, _, _ = specimen_coefficient(xp, outcome, proximity, adjusted=True)
            if np.isfinite(beta):
                betas.append(beta)
        null[b] = np.median(betas)
    p_two = (1 + np.sum(np.abs(null) >= abs(observed_median))) / (N_PERM + 1)
    summary = pd.DataFrame([{
        "outcome": outcome, "proximity": proximity, "n_tumour_specimens": len(obs),
        "observed_median_adjusted_beta": observed_median,
        "n_permutations": N_PERM, "two_sided_fov_permutation_p": p_two,
        "null_q025": np.quantile(null, .025), "null_q50": np.quantile(null, .5), "null_q975": np.quantile(null, .975),
    }])
    return obs, pd.DataFrame({"median_beta": null}), summary


def leave_one_donor_out(coefs: pd.DataFrame) -> pd.DataFrame:
    x = coefs[
        coefs["Sample_type"].eq("T") & coefs["score"].eq("core13")
        & coefs["proximity"].eq("epi_k20") & coefs["model"].eq("adjusted")
    ].dropna(subset=["standardized_beta"])
    rows = []
    for donor in x["Donor"]:
        vals = x.loc[x["Donor"].ne(donor), "standardized_beta"]
        rows.append({"omitted_donor": donor, "n_remaining": len(vals), "median_beta": vals.median(), "n_positive": int((vals > 0).sum())})
    return pd.DataFrame(rows)


def make_figure(d, coefs, tests, perm_null, perm_summary):
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))

    # FOV-binned descriptive relationship; inference remains specimen-level.
    tumour = d[d["Sample_type"].eq("T")].copy()
    tumour["proximity_quintile"] = tumour.groupby(["Slide", "fov"])["epi_k20_fraction"].transform(
        lambda x: pd.qcut(x.rank(method="first"), 5, labels=False, duplicates="drop")
    )
    bins = tumour.groupby(["Sample", "proximity_quintile"], dropna=True).agg(
        proximity=("epi_k20_fraction", "mean"), score=("score_core13", "mean")
    ).reset_index()
    sns.lineplot(data=bins, x="proximity", y="score", units="Sample", estimator=None, color="#B7B7B7", alpha=.55, ax=axes[0, 0])
    sns.regplot(data=bins, x="proximity", y="score", scatter=False, color="#E15759", ax=axes[0, 0])
    axes[0, 0].set_title("a  Tumour specimens: descriptive gradient")
    axes[0, 0].set_xlabel("Local epithelial fraction (20-NN)")
    axes[0, 0].set_ylabel("13-gene cell score")

    main = coefs[(coefs["score"].eq("core13")) & (coefs["proximity"].eq("epi_k20"))]
    sns.stripplot(data=main, x="model", y="standardized_beta", hue="Sample_type", dodge=True, palette={"N": "#4C78A8", "T": "#E15759"}, size=6, ax=axes[0, 1])
    axes[0, 1].axhline(0, color="black", lw=.8, ls="--")
    axes[0, 1].set_title("b  Specimen-level epithelial-gradient slopes")
    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel("Standardized beta")
    axes[0, 1].legend(title="Tissue", frameon=False)

    sens = coefs[(coefs["Sample_type"].eq("T")) & coefs["model"].eq("adjusted")].groupby(["score", "proximity"])["standardized_beta"].median().reset_index()
    sns.barplot(data=sens, x="standardized_beta", y="proximity", hue="score", ax=axes[1, 0])
    axes[1, 0].axvline(0, color="black", lw=.8)
    axes[1, 0].set_title("c  Adjusted score/exposure sensitivity")
    axes[1, 0].set_xlabel("Median tumour-specimen beta")
    axes[1, 0].set_ylabel("")
    axes[1, 0].legend(frameon=False, fontsize=8)

    axes[1, 1].hist(perm_null["median_beta"], bins=38, color="#BAB0AC", edgecolor="white")
    obs = perm_summary.loc[0, "observed_median_adjusted_beta"]
    p = perm_summary.loc[0, "two_sided_fov_permutation_p"]
    axes[1, 1].axvline(obs, color="#E15759", lw=2)
    axes[1, 1].set_title("d  FOV-stratified permutation")
    axes[1, 1].set_xlabel("Median adjusted beta across tumours")
    axes[1, 1].text(.03, .95, f"P={p:.3g}", transform=axes[1, 1].transAxes, va="top")

    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "CommsBio_gate_spatial_gradient.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "CommsBio_gate_spatial_gradient.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    usecols = [
        "Barcode", "Slide", "fov", "Sample", "Donor", "Sample_type", "Quality",
        "Manual_toplevel_pred", "niches", "CenterX_local_px", "CenterY_local_px",
        "CenterX_global_px", "CenterY_global_px", "Mean.PanCK", "Mean.CD45", "Area",
    ]
    meta = pd.read_csv(META, usecols=usecols, low_memory=False)
    meta = meta[meta["Manual_toplevel_pred"].notna() & meta["Manual_toplevel_pred"].ne("QC_fail")].copy()
    schwann = pd.read_csv(SCHWANN, low_memory=False)
    exposures = add_spatial_exposures(meta, set(schwann["Barcode"]))
    extra = meta[["Barcode", "Mean.PanCK", "Mean.CD45", "Area"]].drop_duplicates("Barcode")
    d = schwann.merge(exposures, on=["Barcode", "Slide", "fov"], how="inner").merge(extra, on="Barcode", how="left")
    d["score_core13"] = score_available(d, CORE)
    d["score_noncytokine8"] = score_available(d, NO_CYTOKINES)
    d["score_three_gene"] = score_available(d, THREE_GENE)
    d["myeloid_expr"] = score_available(d, MYELOID)
    d["fibroblast_expr"] = score_available(d, FIBROBLAST)
    d["epithelial_expr"] = score_available(d, EPITHELIAL)
    d["log_library"] = np.log1p(d["nCount_RNA"])
    d["log_panck"] = np.log1p(d["Mean.PanCK"].clip(lower=0))
    d["negative_log_nearest_epi"] = -np.log1p(d["nearest_epi_distance_px"])
    d["negative_log_nearest_epithelial_mass"] = -np.log1p(d["nearest_epithelial_mass_distance_px"])
    d.to_csv(OUT / "schwann_cell_spatial_exposures.csv.gz", index=False)

    coefs = coefficient_table(d)
    tests = specimen_level_tests(coefs)
    perm_obs, perm_null, perm_summary = permute_global_tumour(d)
    loo = leave_one_donor_out(coefs)
    coefs.to_csv(OUT / "specimen_spatial_coefficients.csv", index=False)
    tests.to_csv(OUT / "specimen_level_tests.csv", index=False)
    perm_obs.to_csv(OUT / "permutation_observed_specimens.csv", index=False)
    perm_null.to_csv(OUT / "permutation_null.csv", index=False)
    perm_summary.to_csv(OUT / "permutation_summary.csv", index=False)
    loo.to_csv(OUT / "leave_one_donor_out.csv", index=False)
    make_figure(d, coefs, tests, perm_null, perm_summary)

    focus = tests[(tests["score"].eq("core13")) & (tests["model"].eq("adjusted"))]
    print(f"Spatial exposures calculated for {len(d):,} Schwann cells across {d['Sample'].nunique()} specimens")
    print("\nAdjusted core13 specimen-level tests")
    print(focus.to_string(index=False))
    print("\nFOV-stratified permutation")
    print(perm_summary.to_string(index=False))
    print("\nLeave-one-donor-out range", loo["median_beta"].min(), loo["median_beta"].max())


if __name__ == "__main__":
    main()

