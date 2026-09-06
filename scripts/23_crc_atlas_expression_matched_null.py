#!/usr/bin/env python3
"""Expression-matched three-gene null test in CRC Atlas Schwann cells.

The observed TIMP1-IFITM1-CHI3L1 module was nominated outside CRC Atlas.
This script asks whether its paired tumour-normal effect exceeds effects from
random three-gene sets matched on adjacent-normal Schwann-cell expression,
detection fraction and variance. Patients, not cells, are the inferential unit.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import re
import time

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
H5AD_URL = "https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad"
API_URL = "https://crc.icbi.at/atlas/api/v0.2/data/var"
META_PATH = ROOT / "data/processed/CRC_atlas_schwann_metadata.csv.gz"
CORE_PATH = ROOT / "data/processed/CRC_atlas_schwann_core13_expression.csv.gz"
VAR_STATS_PATH = ROOT / "data/processed/CRC_atlas_var_gene_stats.csv.gz"
CACHE_DIR = ROOT / "data/raw/CRC_atlas_api_matched_pool"
POOL_EXPR_PATH = ROOT / "data/processed/CRC_atlas_schwann_matched_pool_expression.csv.gz"
OUT_DIR = ROOT / "results/tables/fig_matched_null"
FIG_DIR = ROOT / "results/figures"
SOURCE_DIR = ROOT / "results/source_data"

TARGETS = ["TIMP1", "IFITM1", "CHI3L1"]
CORE13 = {
    "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
    "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
}
EXCLUDED_OVERLAP = {"Lee_2020_Nat_Genet"}
N_ATLAS_CELLS = 4_264_929
PRESELECT_PER_TARGET = 35
MATCH_POOL_PER_TARGET = 15
N_RANDOM_SETS = 10_000
SEED = 20260906
MIN_CELLS_PER_SAMPLE = 5


def load_range_reader():
    path = ROOT / "scripts/07_crc_atlas_remote_inspect.py"
    spec = spec_from_file_location("crc_remote", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.HTTPRangeReader


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray([
        x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x)
        for x in values
    ])


def load_or_extract_var_stats() -> pd.DataFrame:
    if VAR_STATS_PATH.exists():
        return pd.read_csv(VAR_STATS_PATH)

    reader = load_range_reader()(H5AD_URL, block_size=4 * 1024 * 1024, max_blocks=12)
    with h5py.File(reader, "r") as h5:
        var = h5["var"]
        frame = pd.DataFrame({
            "gene": decode_strings(var["var_names"][:]),
            "var_index": np.arange(len(var["var_names"]), dtype=np.int32),
            "n_cells_global": var["n_cells"][:],
            "mean_counts_global": var["mean_counts"][:],
            "total_counts_global": var["total_counts"][:],
        })
    frame.to_csv(VAR_STATS_PATH, index=False, compression="gzip")
    return frame


def robust_z(values: pd.Series) -> pd.Series:
    median = values.median()
    mad = np.median(np.abs(values - median))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale == 0:
        scale = values.std(ddof=0)
    return (values - median) / scale


def preselect_candidates(var_stats: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    d = var_stats.copy()
    d = d[d["gene"].str.match(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", na=False)]
    d = d[~d["gene"].duplicated(keep=False)]
    d = d[d["n_cells_global"].between(50, N_ATLAS_CELLS - 1)]
    d = d[~d["gene"].isin(CORE13)]
    d = d[~d["gene"].str.upper().str.match(r"^(MT-|RPL|RPS)")]

    d["global_log_mean"] = np.log1p(d["mean_counts_global"].clip(lower=0))
    p = (d["n_cells_global"] + 0.5) / (N_ATLAS_CELLS + 1.0)
    d["global_logit_detection"] = np.log(p / (1 - p))

    target_rows = var_stats.set_index("gene").loc[TARGETS].copy()
    target_rows["global_log_mean"] = np.log1p(target_rows["mean_counts_global"].clip(lower=0))
    tp = (target_rows["n_cells_global"] + 0.5) / (N_ATLAS_CELLS + 1.0)
    target_rows["global_logit_detection"] = np.log(tp / (1 - tp))

    combined_mean = pd.concat([d["global_log_mean"], target_rows["global_log_mean"]])
    combined_det = pd.concat([d["global_logit_detection"], target_rows["global_logit_detection"]])
    mean_median, mean_scale = combined_mean.median(), 1.4826 * np.median(np.abs(combined_mean - combined_mean.median()))
    det_median, det_scale = combined_det.median(), 1.4826 * np.median(np.abs(combined_det - combined_det.median()))
    d["z_global_mean"] = (d["global_log_mean"] - mean_median) / mean_scale
    d["z_global_detection"] = (d["global_logit_detection"] - det_median) / det_scale
    target_rows["z_global_mean"] = (target_rows["global_log_mean"] - mean_median) / mean_scale
    target_rows["z_global_detection"] = (target_rows["global_logit_detection"] - det_median) / det_scale

    rows = []
    selected = set()
    for target in TARGETS:
        dist = np.sqrt(
            (d["z_global_mean"] - target_rows.loc[target, "z_global_mean"]) ** 2
            + (d["z_global_detection"] - target_rows.loc[target, "z_global_detection"]) ** 2
        )
        nearest = d.assign(global_preselection_distance=dist).nsmallest(PRESELECT_PER_TARGET, "global_preselection_distance")
        nearest = nearest.assign(target=target, preselection_rank=np.arange(1, len(nearest) + 1))
        rows.append(nearest)
        selected.update(nearest["gene"])
    out = pd.concat(rows, ignore_index=True)
    return out, sorted(selected)


def load_decoder():
    path = ROOT / "scripts/19_crc_atlas_schwann_expression.py"
    spec = spec_from_file_location("crc_expr", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.decode_single_float_column


def safe_cache_name(var_index: int, gene: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]", "_", gene)
    return f"{var_index:05d}_{token}.fbs"


def fetch_one_gene(gene: str, var_index: int, decoder, rows: np.ndarray) -> tuple[str, np.ndarray]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / safe_cache_name(var_index, gene)
    payload = None
    if path.exists() and path.stat().st_size > 1_000_000:
        payload = path.read_bytes()
    else:
        last_error = None
        for attempt in range(5):
            try:
                response = requests.get(
                    API_URL,
                    params={"var:var_names": gene},
                    headers={
                        "Accept": "application/octet-stream",
                        "Accept-Encoding": "gzip, deflate",
                        "User-Agent": "Mozilla/5.0",
                    },
                    timeout=120,
                )
                response.raise_for_status()
                if len(response.content) < 1_000_000:
                    raise RuntimeError(f"short response ({len(response.content)} bytes)")
                payload = response.content
                tmp = path.with_suffix(path.suffix + ".part")
                tmp.write_bytes(payload)
                tmp.replace(path)
                break
            except Exception as exc:  # network retries are bounded and reported
                last_error = exc
                time.sleep(1.5 * (attempt + 1))
        if payload is None:
            raise RuntimeError(f"Failed to fetch {gene}: {last_error}")

    values, returned_index = decoder(payload)
    if len(values) != N_ATLAS_CELLS:
        raise RuntimeError(f"Atlas row count changed for {gene}: {len(values)}")
    if returned_index != var_index:
        raise RuntimeError(f"Var-index mismatch for {gene}: expected {var_index}, received {returned_index}")
    return gene, np.asarray(values[rows], dtype=np.float32)


def load_or_fetch_pool_expression(candidates: list[str], var_stats: pd.DataFrame) -> pd.DataFrame:
    metadata = pd.read_csv(META_PATH)
    rows = metadata["obs_index"].to_numpy(np.int64)
    existing = pd.DataFrame(index=np.arange(len(rows)))
    if POOL_EXPR_PATH.exists():
        existing = pd.read_csv(POOL_EXPR_PATH)

    missing = [g for g in candidates if g not in existing.columns]
    if missing:
        lookup = var_stats.set_index("gene")["var_index"].astype(int).to_dict()
        decoder = load_decoder()
        fetched = {}
        # The public server is reliable with compressed responses but can close
        # concurrent uncompressed streams; two workers keep load conservative.
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {
                executor.submit(fetch_one_gene, gene, lookup[gene], decoder, rows): gene
                for gene in missing
            }
            for i, future in enumerate(as_completed(futures), start=1):
                gene, values = future.result()
                fetched[gene] = values
                print(f"Fetched {i}/{len(missing)}: {gene}", flush=True)
        new = pd.DataFrame(fetched)
        existing = pd.concat([existing.reset_index(drop=True), new.reset_index(drop=True)], axis=1)
        existing.to_csv(POOL_EXPR_PATH, index=False, compression="gzip")
    return existing[candidates]


def normal_feature_table(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    normal_mask = metadata["sample_type"].eq("adjacent normal").to_numpy()
    x = expression.loc[normal_mask]
    return pd.DataFrame({
        "gene": x.columns,
        "normal_mean": x.mean(axis=0).to_numpy(),
        "normal_detection": x.gt(0).mean(axis=0).to_numpy(),
        "normal_variance": x.var(axis=0, ddof=1).to_numpy(),
    })


def exact_match_pools(
    features: pd.DataFrame,
    pool_size: int = MATCH_POOL_PER_TARGET,
    excluded_gene_pattern: str | None = None,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    d = features.copy()
    d["t_mean"] = np.log1p(d["normal_mean"].clip(lower=0))
    p = d["normal_detection"].clip(1e-5, 1 - 1e-5)
    d["t_detection"] = np.log(p / (1 - p))
    d["t_variance"] = np.log1p(d["normal_variance"].clip(lower=0))
    for col in ["t_mean", "t_detection", "t_variance"]:
        d[f"z_{col}"] = robust_z(d[col])

    indexed = d.set_index("gene")
    rows = []
    pools = {}
    for target in TARGETS:
        available = d[~d["gene"].isin(TARGETS)].copy()
        if excluded_gene_pattern is not None:
            available = available[
                ~available["gene"].str.match(excluded_gene_pattern, case=False, na=False)
            ]
        available["match_distance"] = np.sqrt(sum(
            (available[f"z_{col}"] - indexed.loc[target, f"z_{col}"]) ** 2
            for col in ["t_mean", "t_detection", "t_variance"]
        ))
        nearest = available.nsmallest(pool_size, "match_distance").copy()
        nearest["target"] = target
        nearest["match_rank"] = np.arange(1, len(nearest) + 1)
        for col in ["t_mean", "t_detection", "t_variance"]:
            nearest[f"target_z_{col}"] = indexed.loc[target, f"z_{col}"]
            nearest[f"standardized_difference_{col}"] = nearest[f"z_{col}"] - indexed.loc[target, f"z_{col}"]
        rows.append(nearest)
        pools[target] = nearest["gene"].tolist()
    return pd.concat(rows, ignore_index=True), pools


def depth_adjust_all(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    adjusted = pd.DataFrame(index=expression.index, columns=expression.columns, dtype=np.float32)
    for _, idx in metadata.groupby("study_id", observed=True).groups.items():
        pos = np.asarray(list(idx), dtype=int)
        part = metadata.loc[pos]
        design = np.column_stack([
            np.ones(len(part)),
            np.log1p(part["n_counts"].to_numpy(float)),
            np.log1p(part["n_genes"].to_numpy(float)),
        ])
        y = expression.iloc[pos].to_numpy(float)
        coef = np.linalg.lstsq(design, y, rcond=None)[0]
        adjusted.iloc[pos] = (y - design @ coef + y.mean(axis=0)).astype(np.float32)
    return adjusted


def patient_gene_deltas(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    adjusted = depth_adjust_all(expression, metadata)
    keys = ["study_id", "patient_id", "sample_id", "sample_type"]
    sample_scores = pd.concat([metadata[keys], adjusted], axis=1).groupby(keys, observed=True).mean().reset_index()
    sample_n = metadata.groupby(keys, observed=True).size().rename("n_cells").reset_index()
    sample_scores = sample_scores.merge(sample_n, on=keys, how="left")
    sample_scores = sample_scores[sample_scores["n_cells"] >= MIN_CELLS_PER_SAMPLE]
    sample_scores = sample_scores[sample_scores["sample_type"].isin(["primary tumor", "adjacent normal"])]
    patient = sample_scores.groupby(["study_id", "patient_id", "sample_type"], observed=True)[expression.columns.tolist()].mean()

    rows = []
    for gene in expression.columns:
        wide = patient[gene].unstack("sample_type").dropna(subset=["primary tumor", "adjacent normal"])
        delta = wide["primary tumor"] - wide["adjacent normal"]
        rows.append(delta.rename(gene))
    out = pd.concat(rows, axis=1)
    out.index.names = ["study_id", "patient_id"]
    return out.reset_index()


def random_set_test(deltas: pd.DataFrame, pools: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    indexed = deltas.set_index(["study_id", "patient_id"])
    observed_values = indexed[TARGETS].mean(axis=1)
    observed_median = float(observed_values.median())
    observed_positive = int((observed_values > 0).sum())

    rng = np.random.default_rng(SEED)
    rows = []
    for iteration in range(1, N_RANDOM_SETS + 1):
        for _ in range(100):
            chosen = [rng.choice(pools[target]) for target in TARGETS]
            if len(set(chosen)) == len(chosen):
                break
        else:
            raise RuntimeError("Unable to draw three unique matched genes")
        values = indexed[chosen].mean(axis=1)
        rows.append({
            "iteration": iteration,
            "TIMP1_match": chosen[0],
            "IFITM1_match": chosen[1],
            "CHI3L1_match": chosen[2],
            "median_paired_delta": float(values.median()),
            "n_positive_patients": int((values > 0).sum()),
        })
    null = pd.DataFrame(rows)
    p_delta = (1 + (null["median_paired_delta"] >= observed_median).sum()) / (len(null) + 1)
    p_positive = (1 + (null["n_positive_patients"] >= observed_positive).sum()) / (len(null) + 1)
    summary = pd.DataFrame([{
        "seed": SEED,
        "n_random_sets": len(null),
        "n_paired_patients": len(indexed),
        "n_studies": deltas["study_id"].nunique(),
        "observed_genes": "-".join(TARGETS),
        "observed_median_paired_delta": observed_median,
        "observed_n_positive_patients": observed_positive,
        "null_median_delta_q025": null["median_paired_delta"].quantile(0.025),
        "null_median_delta_q50": null["median_paired_delta"].quantile(0.5),
        "null_median_delta_q975": null["median_paired_delta"].quantile(0.975),
        "empirical_p_median_delta_one_sided": p_delta,
        "null_n_positive_q025": null["n_positive_patients"].quantile(0.025),
        "null_n_positive_q50": null["n_positive_patients"].quantile(0.5),
        "null_n_positive_q975": null["n_positive_patients"].quantile(0.975),
        "empirical_p_n_positive_one_sided": p_positive,
        "observed_effect_percentile": 100 * (null["median_paired_delta"] < observed_median).mean(),
    }])
    return null, summary


def make_figure(null: pd.DataFrame, summary: pd.DataFrame, matches: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8))
    s = summary.iloc[0]

    axes[0].hist(null["median_paired_delta"], bins=45, color="#B8B8B8", edgecolor="white")
    axes[0].axvline(s["observed_median_paired_delta"], color="#D1495B", lw=2.2)
    axes[0].set(
        title="a  Matched-set effect null",
        xlabel="Median paired tumour−normal difference",
        ylabel="Random three-gene sets",
    )
    axes[0].text(
        0.03, 0.96,
        f"Observed={s['observed_median_paired_delta']:.3f}\nEmpirical P={s['empirical_p_median_delta_one_sided']:.4f}",
        transform=axes[0].transAxes, va="top",
    )

    bins = np.arange(null["n_positive_patients"].min() - 0.5, null["n_positive_patients"].max() + 1.5)
    axes[1].hist(null["n_positive_patients"], bins=bins, color="#76A5AF", edgecolor="white")
    axes[1].axvline(s["observed_n_positive_patients"], color="#D1495B", lw=2.2)
    axes[1].set(
        title="b  Directional-consistency null",
        xlabel="Patients with tumour > normal",
        ylabel="Random three-gene sets",
    )
    axes[1].text(
        0.03, 0.96,
        f"Observed={int(s['observed_n_positive_patients'])}/{int(s['n_paired_patients'])}\nEmpirical P={s['empirical_p_n_positive_one_sided']:.4f}",
        transform=axes[1].transAxes, va="top",
    )

    feature_labels = {
        "t_mean": "Mean",
        "t_detection": "Detection",
        "t_variance": "Variance",
    }
    plot_rows = []
    for _, row in matches.iterrows():
        for feature, label in feature_labels.items():
            plot_rows.append({
                "target": row["target"],
                "feature": label,
                "absolute_standardized_difference": abs(row[f"standardized_difference_{feature}"]),
            })
    quality = pd.DataFrame(plot_rows)
    sns.boxplot(
        data=quality, x="feature", y="absolute_standardized_difference", hue="target",
        showfliers=False, palette="Set2", ax=axes[2],
    )
    sns.stripplot(
        data=quality, x="feature", y="absolute_standardized_difference", hue="target",
        dodge=True, size=2.5, alpha=0.55, palette="Set2", legend=False, ax=axes[2],
    )
    axes[2].axhline(0.5, color="#666666", ls="--", lw=1)
    axes[2].set(
        title="c  Adjacent-normal match quality",
        xlabel="Matched feature",
        ylabel="Absolute standardized difference",
    )
    axes[2].legend(title="Target", frameon=False, fontsize=7, title_fontsize=8)

    sns.despine(fig)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "FIG_expression_matched_random_set_null.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "FIG_expression_matched_random_set_null.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    var_stats = load_or_extract_var_stats()
    preselection, candidates = preselect_candidates(var_stats)
    preselection.to_csv(OUT_DIR / "candidate_preselection.csv", index=False)

    metadata = pd.read_csv(META_PATH)
    keep = ~metadata["study_id"].isin(EXCLUDED_OVERLAP)
    metadata = metadata.loc[keep].reset_index(drop=True)

    pool_expression_all = load_or_fetch_pool_expression(candidates, var_stats)
    pool_expression = pool_expression_all.loc[keep.to_numpy()].reset_index(drop=True)
    core = pd.read_csv(CORE_PATH, usecols=TARGETS).loc[keep].reset_index(drop=True)
    expression = pd.concat([core, pool_expression], axis=1)

    features = normal_feature_table(expression, metadata)
    matches, pools = exact_match_pools(features)
    features.to_csv(OUT_DIR / "adjacent_normal_gene_features.csv", index=False)
    matches.to_csv(OUT_DIR / "target_matched_gene_pools.csv", index=False)

    # Retain every preselected candidate here so prespecified pool-size and
    # immune-receptor-exclusion sensitivity analyses use the same patient set.
    analysis_expression = expression[TARGETS + candidates]
    deltas = patient_gene_deltas(analysis_expression, metadata)
    deltas.to_csv(OUT_DIR / "paired_patient_gene_deltas.csv.gz", index=False, compression="gzip")
    null, summary = random_set_test(deltas, pools)
    null.to_csv(OUT_DIR / "expression_matched_random_set_null.csv.gz", index=False, compression="gzip")
    summary.to_csv(OUT_DIR / "expression_matched_random_set_test.csv", index=False)
    null.to_csv(SOURCE_DIR / "FigureS2_expression_matched_random_set_null.csv", index=False)
    matches.to_csv(SOURCE_DIR / "FigureS2_matching_diagnostics.csv", index=False)

    sensitivity_rows = []
    sensitivity_specs = [
        ("nearest_5", 5, None),
        ("nearest_10", 10, None),
        ("nearest_15_primary", 15, None),
        ("nearest_15_excluding_IG_TR", 15, r"^(IG[HKL]|TR[ABDG])"),
    ]
    for specification, pool_size, excluded_pattern in sensitivity_specs:
        _, sensitivity_pools = exact_match_pools(
            features,
            pool_size=pool_size,
            excluded_gene_pattern=excluded_pattern,
        )
        _, sensitivity_summary = random_set_test(deltas, sensitivity_pools)
        sensitivity_summary.insert(0, "specification", specification)
        sensitivity_summary.insert(1, "matched_candidates_per_target", pool_size)
        sensitivity_rows.append(sensitivity_summary)
    sensitivity = pd.concat(sensitivity_rows, ignore_index=True)
    sensitivity.to_csv(OUT_DIR / "expression_matched_random_set_sensitivity.csv", index=False)
    sensitivity.to_csv(SOURCE_DIR / "FigureS2_matched_null_sensitivity.csv", index=False)
    make_figure(null, summary, matches)

    print("\nPrimary result")
    print(summary.to_string(index=False))
    print("\nMatched pools")
    print(matches.groupby("target").agg(
        n_candidates=("gene", "size"),
        median_distance=("match_distance", "median"),
        max_distance=("match_distance", "max"),
    ).to_string())


if __name__ == "__main__":
    main()
