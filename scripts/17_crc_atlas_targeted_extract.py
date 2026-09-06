#!/usr/bin/env python3
"""Extract only CRC-Atlas Schwann metadata through remote HDF5 byte ranges.

The 32.7-GB source file is never downloaded.  Categorical observation columns
are represented by small code vectors and dictionaries, so they can be sliced
directly and saved for subsequent expression requests through the cellxgene
API.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from importlib.util import module_from_spec, spec_from_file_location


ROOT = Path(__file__).resolve().parents[1]
URL = "https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad"
OUT = ROOT / "data/processed/CRC_atlas_schwann_metadata.csv.gz"
FIELDS = [
    "dataset",
    "medical_condition",
    "cancer_type",
    "sample_id",
    "sample_type",
    "tumor_source",
    "sample_tissue",
    "anatomic_region",
    "patient_id",
    "study_id",
    "atlas_cell_type_coarse",
    "atlas_cell_type_middle",
    "n_counts",
    "n_genes",
]


def load_range_reader():
    path = ROOT / "scripts/07_crc_atlas_remote_inspect.py"
    spec = spec_from_file_location("crc_remote", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.HTTPRangeReader


def decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O"}:
        return np.asarray([
            x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else x
            for x in values
        ])
    return values


def read_column(obs: h5py.Group, name: str, rows: np.ndarray | None = None) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:] if rows is None else obj["codes"][rows]
        categories = decode_array(obj["categories"][:])
        result = np.empty(len(codes), dtype=object)
        result[:] = None
        valid = codes >= 0
        result[valid] = categories[codes[valid]]
        return result
    values = obj[:] if rows is None else obj[rows]
    return decode_array(values)


def main() -> None:
    HTTPRangeReader = load_range_reader()
    remote = HTTPRangeReader(URL, block_size=4 * 1024 * 1024, max_blocks=12)
    with h5py.File(remote, "r") as h5:
        obs = h5["obs"]
        middle = obs["atlas_cell_type_middle"]
        categories = decode_array(middle["categories"][:])
        matches = np.flatnonzero(categories == "Schwann cell")
        if len(matches) != 1:
            raise RuntimeError(f"Expected one Schwann category; found {matches.tolist()}")
        target_code = int(matches[0])
        codes = middle["codes"][:]
        rows = np.flatnonzero(codes == target_code).astype(np.int64)
        print(f"Found {len(rows):,} Schwann cells among {len(codes):,} atlas cells")

        frame = pd.DataFrame({"obs_index": rows})
        for field in FIELDS:
            frame[field] = read_column(obs, field, rows)

        OUT.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(OUT, index=False, compression="gzip")
        summary = {
            "source": URL,
            "n_atlas_cells": int(len(codes)),
            "n_schwann_cells": int(len(rows)),
            "n_patients": int(frame["patient_id"].nunique()),
            "n_samples": int(frame["sample_id"].nunique()),
            "n_studies": int(frame["study_id"].nunique()),
            "conditions": frame["medical_condition"].value_counts(dropna=False).to_dict(),
        }
        OUT.with_suffix("").with_suffix(".json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
