from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_genes(path: Path) -> tuple[dict[str, list[str]], list[str]]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    return cfg, sorted({gene for genes in cfg.values() for gene in genes})


def extract_rows(matrix_path: Path, wanted: set[str]) -> tuple[pd.DataFrame, pd.Series]:
    opener = gzip.open if matrix_path.suffix == ".gz" else open
    with opener(matrix_path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = [x.strip('"') for x in header]
        rows: dict[str, np.ndarray] = {}
        library: np.ndarray | None = None
        for line in handle:
            gene_text, values_text = line.rstrip("\n").split("\t", 1)
            gene = gene_text.strip('"')
            values = np.fromstring(values_text, sep="\t", dtype=np.float32)
            if library is None:
                library = np.zeros_like(values, dtype=np.float64)
            library += values
            if gene in wanted:
                rows[gene] = values
    if not rows:
        raise ValueError(f"No requested genes found in {matrix_path}")
    expected = len(next(iter(rows.values())))
    if len(cells) == expected + 1 and cells[0] in {"", "Index", "Gene", "gene"}:
        cells = cells[1:]
    if len(cells) != expected:
        raise ValueError(
            f"Header has {len(cells)} fields but rows have {expected} values"
        )
    if library is None:
        raise ValueError(f"Empty matrix: {matrix_path}")
    return pd.DataFrame(rows, index=cells), pd.Series(library, index=cells, name="library")


def summarize(dataset: str, expression: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    annotation = annotation.copy()
    annotation.index = annotation["Index"].astype(str)
    common = expression.index.intersection(annotation.index)
    expression = expression.loc[common]
    annotation = annotation.loc[common]
    joined = annotation.join(expression)
    group_cols = ["Patient", "Class", "Cell_type", "Cell_subtype"]
    value_cols = list(expression.columns)
    means = joined.groupby(group_cols, observed=True)[value_cols].mean()
    detected = joined.assign(**{g: joined[g].gt(0) for g in value_cols}).groupby(
        group_cols, observed=True
    )[value_cols].mean()
    counts = joined.groupby(group_cols, observed=True).size().rename("n_cells")
    out = means.add_suffix("__mean").join(detected.add_suffix("__pct_detected")).join(counts)
    out.insert(0, "dataset", dataset)
    return out.reset_index()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--genesets", type=Path, default=Path("config/genesets.yaml"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    gene_sets, genes = load_genes(args.genesets)
    raw_expression, library = extract_rows(args.matrix, set(genes))
    expression = np.log1p(raw_expression.divide(library.clip(lower=1), axis=0) * 10_000)
    annotation = pd.read_csv(args.annotation, sep="\t")
    summary = summarize(args.dataset, expression, annotation)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, index=False)

    common = expression.index.intersection(annotation["Index"].astype(str))
    annotation_by_cell = annotation.set_index(annotation["Index"].astype(str)).loc[common]
    reactive_genes = [g for g in gene_sets["reactive_egc_conserved"] if g in expression]
    reactive_score = expression.loc[common, reactive_genes].mean(axis=1)
    score_data = annotation_by_cell[["Patient", "Class", "Cell_type", "Cell_subtype"]].copy()
    score_data["reactive_egc_score"] = reactive_score
    score_summary = score_data.groupby(
        ["Patient", "Class", "Cell_type", "Cell_subtype"], observed=True
    )["reactive_egc_score"].agg(["mean", "median", "count"]).reset_index()
    score_summary.to_csv(args.output.with_name(args.output.stem + "_reactive_scores.csv"), index=False)

    report = {
        "dataset": args.dataset,
        "n_annotation_cells": int(len(annotation)),
        "n_matched_cells": int(expression.index.intersection(annotation["Index"]).size),
        "requested_genes": genes,
        "found_genes": sorted(expression.columns),
        "missing_genes": sorted(set(genes) - set(expression.columns)),
    }
    report_path = args.output.with_suffix(".json")
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
