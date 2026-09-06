#!/usr/bin/env python3
"""Create the supplementary study-selection figure and source-data manifest.

Figure contract
---------------
Core conclusion: the analysis separates programme derivation, independent
population validation, orthogonal spatial-platform validation, and supportive
cohorts while keeping patients, donors, or samples as biological replicates.
Archetype: schematic-led composite.
Target: Communications Biology supplementary figure, 183 mm wide.
Backend: Python/matplotlib only.
Source data: every displayed count is written to a CSV beside the figure data.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "manuscript_figures"
DATA_DIR = ROOT / "results" / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

INK = "#243447"
BLUE = "#4C78A8"
BLUE_LIGHT = "#E8F0F7"
TEAL = "#2A9D8F"
TEAL_LIGHT = "#E5F4F1"
GOLD = "#C58A24"
GOLD_LIGHT = "#FBF1DD"
GREY = "#6B7280"
GREY_LIGHT = "#F2F4F7"
RED = "#B24A3B"
RED_LIGHT = "#F8EAE7"


SELECTION_ROWS = [
    {
        "stage": "Programme derivation",
        "resource": "GSE231802",
        "input": "Tumour-conditioned mouse enteric glia; 3 time points; n=4 per group/time",
        "filter": "log2 fold change >0.5 and BH FDR <0.05 at >=2 time points",
        "output": "33-gene programme; 13 genes shared by both CosMx panels",
        "biological_unit": "Culture replicate",
        "role": "Programme generator",
    },
    {
        "stage": "Primary validation",
        "resource": "CRC Atlas",
        "input": "4,264,929 cells; 15,352 annotated Schwann cells; 279 patients; 560 samples; 32 studies",
        "filter": "Exclude Lee_2020_Nat_Genet; primary tumour and adjacent normal; >=5 Schwann cells per tissue",
        "output": "44 paired patients from 9 studies",
        "biological_unit": "Paired patient",
        "role": "Independent population validation",
    },
    {
        "stage": "Orthogonal validation",
        "resource": "GSE303070",
        "input": "846,469 cells across 352 fields of view",
        "filter": "Author QC labels; author Schwann annotation; specimen pseudobulk",
        "output": "803,109 QC-passing cells; 29,829 Schwann cells; 24 specimens; 16 donors; 8 paired donors",
        "biological_unit": "Paired donor",
        "role": "CosMx validation and module nomination",
    },
    {
        "stage": "Supportive validation",
        "resource": "GSE132465",
        "input": "305 annotated enteric-glial cells",
        "filter": ">=3 cells per patient-tissue summary",
        "output": "12 summaries: 3 tumour and 9 normal",
        "biological_unit": "Patient-tissue summary",
        "role": "Supportive single-cell cohort",
    },
    {
        "stage": "Supportive validation",
        "resource": "GSE144735",
        "input": "463 annotated enteric-glial cells",
        "filter": ">=3 cells per patient-tissue summary",
        "output": "9 summaries: 3 tumour and 6 normal",
        "biological_unit": "Patient-tissue summary",
        "role": "Supportive single-cell cohort",
    },
]


def box(ax, xy, wh, title, body, edge, face, title_color=None, lw=1.0):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + 0.018, y + h - 0.025, title, ha="left", va="top", fontsize=7.4,
            fontweight="bold", color=title_color or edge)
    ax.text(x + 0.018, y + h - 0.075, body, ha="left", va="top", fontsize=6.0,
            color=INK, linespacing=1.25)


def arrow(ax, start, end, color=GREY):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.9,
            color=color,
            shrinkA=2,
            shrinkB=2,
        )
    )


def write_source_data():
    out = DATA_DIR / "FigureS1_study_selection.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SELECTION_ROWS[0]))
        writer.writeheader()
        writer.writerows(SELECTION_ROWS)
    return out


def build_figure():
    width_mm = 183
    height_mm = 118
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.0, 0.985, "a", fontsize=8, fontweight="bold", va="top", color=INK)
    ax.text(0.035, 0.985, "Programme derivation and prespecification", fontsize=7.5,
            fontweight="bold", va="top", color=INK)

    box(ax, (0.035, 0.775), (0.25, 0.13), "GSE231802 perturbation",
        "Mouse enteric glia\n3 conditions x 3 time points\nn = 4 per group/time",
        BLUE, BLUE_LIGHT)
    box(ax, (0.375, 0.775), (0.25, 0.13), "Time-consistent response",
        "log2FC > 0.5 and BH FDR < 0.05\nat two or more time points\n33 genes retained",
        BLUE, BLUE_LIGHT)
    box(ax, (0.715, 0.775), (0.25, 0.13), "Frozen cross-platform scaffold",
        "13 genes measured by both\nCosMx panel versions\nFrozen before human tests",
        TEAL, TEAL_LIGHT)
    arrow(ax, (0.285, 0.84), (0.375, 0.84), BLUE)
    arrow(ax, (0.625, 0.84), (0.715, 0.84), TEAL)

    ax.text(0.0, 0.69, "b", fontsize=8, fontweight="bold", va="top", color=INK)
    ax.text(0.035, 0.69, "Human validation resources and analysis units", fontsize=7.5,
            fontweight="bold", va="top", color=INK)

    box(ax, (0.035, 0.49), (0.29, 0.145), "CRC Atlas input",
        "4,264,929 cells\n15,352 annotated Schwann cells\n279 patients, 560 samples, 32 studies",
        GOLD, GOLD_LIGHT)
    box(ax, (0.355, 0.49), (0.29, 0.145), "Eligibility and independence",
        "Excluded Lee_2020_Nat_Genet\nPrimary tumour + adjacent normal\n>=5 Schwann cells in each tissue",
        RED, RED_LIGHT)
    box(ax, (0.675, 0.49), (0.29, 0.145), "Primary validation set",
        "44 paired patients\n9 independent studies\nPatient pair = inferential unit",
        GOLD, GOLD_LIGHT)
    arrow(ax, (0.325, 0.562), (0.355, 0.562), GOLD)
    arrow(ax, (0.645, 0.562), (0.675, 0.562), GOLD)

    box(ax, (0.035, 0.255), (0.29, 0.16), "CosMx GSE303070",
        "846,469 cells; 352 fields of view\n803,109 cells after author QC\n29,829 Schwann cells; 24 specimens\n16 donors; 8 paired donors",
        TEAL, TEAL_LIGHT)
    box(ax, (0.355, 0.255), (0.29, 0.16), "GSE132465",
        "305 enteric-glial cells\n>=3 cells per patient-tissue summary\n12 summaries: 3 tumour, 9 normal",
        BLUE, BLUE_LIGHT)
    box(ax, (0.675, 0.255), (0.29, 0.16), "GSE144735",
        "463 enteric-glial cells\n>=3 cells per patient-tissue summary\n9 summaries: 3 tumour, 6 normal",
        BLUE, BLUE_LIGHT)

    box(ax, (0.035, 0.045), (0.93, 0.135), "Inference hierarchy",
        "Programme generator: culture replicate   |   Primary test: paired patient   |   "
        "Orthogonal test: paired donor   |   Supportive tests: patient-tissue summary\n"
        "Cells contribute to aggregate scores; they do not enter inferential tests as independent replicates.",
        GREY, GREY_LIGHT, title_color=INK)
    arrow(ax, (0.18, 0.255), (0.25, 0.18), TEAL)
    arrow(ax, (0.50, 0.255), (0.50, 0.18), BLUE)
    arrow(ax, (0.82, 0.255), (0.75, 0.18), BLUE)

    fig.subplots_adjust(left=0.025, right=0.985, top=0.985, bottom=0.02)
    stem = FIG_DIR / "FigureS1_study_selection"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return stem


if __name__ == "__main__":
    print(write_source_data())
    print(build_figure())
