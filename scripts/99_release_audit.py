#!/usr/bin/env python3
"""Fail when a release contains missing evidence, secrets, local paths, or bulky raw files."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {".git", ".mplconfig", ".uv-cache", ".uv-python", ".venv", "__pycache__", "logs"}
TEXT_SUFFIXES = {".cff", ".csv", ".json", ".md", ".py", ".r", ".toml", ".tsv", ".txt", ".yaml", ".yml"}
BANNED_SUFFIXES = {".fbs", ".h5", ".h5ad", ".h5seurat", ".loom", ".mtx", ".parquet", ".rds", ".zip"}
PATTERNS = {
    "home path": re.compile(r"/(Users|home)/[^/\s]+/"),
    "GitHub token": re.compile(r"\b(?:ghp|github_pat)_[A-Za-z0-9_]+"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "authorization header": re.compile(r"Authorization\s*:\s*(?:Bearer|token)\s+\S+", re.I),
    "unresolved repository placeholder": re.compile(
        "REPOSITORY" + "_OWNER|AUTHOR INPUT" + " REQUIRED", re.I
    ),
}

REQUIRED_SOURCE = {
    "Figure1_reactive_program.csv",
    "Figure2_CRC_atlas_paired_patients.csv",
    "Figure2_CRC_atlas_studies.csv",
    "Figure2_CRC_atlas_threshold_sensitivity.csv",
    "Figure2_CRC_atlas_leave_one_study_out.csv",
    "Figure3_paired_scores.csv",
    "Figure3_paired_gene_deltas.csv",
    "Figure3_gene_tests.csv",
    "Figure4_GSE132465_scores.csv",
    "Figure4_GSE144735_scores.csv",
    "Figure4_scrna_threshold_sensitivity.csv",
    "Figure4_spatial_boundary.csv",
    "FigureS1_study_selection.csv",
    "FigureS2_expression_matched_random_set_null.csv",
    "FigureS2_matching_diagnostics.csv",
    "FigureS2_matched_null_sensitivity.csv",
    "source_data_manifest.csv",
}


def main() -> None:
    problems: list[str] = []
    gene_sets: dict[str, list[str]] = {}
    current: str | None = None
    for raw in (ROOT / "config/genesets.yaml").read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line:
            continue
        if not line.startswith(" ") and line.endswith(":"):
            current = line[:-1]
            gene_sets[current] = []
        elif current and line.lstrip().startswith("- "):
            gene_sets[current].append(line.lstrip()[2:].strip())
    expected = {
        "perturbation_derived_program_33": 33,
        "cross_platform_scaffold_13": 13,
        "convergent_human_module_3": 3,
    }
    for name, count in expected.items():
        values = gene_sets.get(name, [])
        if len(values) != count or len(set(values)) != count:
            problems.append(f"{name}: expected {count} unique genes, found {len(set(values))}")

    source_dir = ROOT / "results/source_data"
    missing = sorted(REQUIRED_SOURCE - {p.name for p in source_dir.glob("*.csv")})
    if missing:
        problems.append("missing source-data files: " + ", ".join(missing))

    for path in ROOT.rglob("*"):
        if not path.is_file() or IGNORED_DIRS.intersection(path.relative_to(ROOT).parts):
            continue
        rel = path.relative_to(ROOT)
        if path.suffix.lower() in BANNED_SUFFIXES:
            problems.append(f"excluded raw/cache format present: {rel}")
        if path.stat().st_size > 50 * 1024 * 1024:
            problems.append(f"file exceeds 50 MiB: {rel}")
        if path.suffix.lower() in TEXT_SUFFIXES or path.name in {"LICENSE", "Makefile"}:
            content = path.read_text(encoding="utf-8", errors="ignore")
            for label, pattern in PATTERNS.items():
                if pattern.search(content):
                    problems.append(f"{label} in {rel}")

    if problems:
        raise SystemExit("Release audit failed:\n- " + "\n- ".join(sorted(set(problems))))
    print("Release audit passed")
    print("Fixed gene sets: 33, 13, and 3 unique genes")
    print(f"Source-data files: {len(REQUIRED_SOURCE)} required files present")


if __name__ == "__main__":
    main()
