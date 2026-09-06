#!/usr/bin/env python3
"""Download 13 CRC-Atlas expression columns and retain Schwann cells only.

CELLxGENE serves each gene as a FlatBuffer matrix.  The response contains one
float32 vector over all atlas observations; only pre-identified Schwann indices
are retained.  Raw responses are cached so reruns do not redownload data.
"""

from __future__ import annotations

from pathlib import Path
import struct

import flatbuffers
from flatbuffers import encode, number_types, packer
from flatbuffers.table import Table
import numpy as np
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data/processed/CRC_atlas_schwann_metadata.csv.gz"
CACHE = ROOT / "data/raw/CRC_atlas_api"
OUT = ROOT / "data/processed/CRC_atlas_schwann_core13_expression.csv.gz"
URL = "https://crc.icbi.at/atlas/api/v0.2/data/var"
CORE = [
    "MT2A", "SOD2", "CCL2", "IER3", "CXCL10", "IL6", "LCN2", "SLPI",
    "CHI3L1", "CSF3", "CCL11", "TIMP1", "IFITM1",
]
N_ATLAS_CELLS = 4_264_929


def table_scalar(tab: Table, slot: int, flag):
    offset = number_types.UOffsetTFlags.py_type(tab.Offset(slot))
    return tab.Get(flag, offset + tab.Pos) if offset else 0


def union_vector(tab: Table, union_slot: int, type_flag):
    offset = number_types.UOffsetTFlags.py_type(tab.Offset(union_slot))
    union = Table(bytearray(), 0)
    tab.Union(union, offset)
    data_offset = number_types.UOffsetTFlags.py_type(union.Offset(4))
    return union.GetVectorAsNumpy(type_flag, data_offset)


def decode_single_float_column(payload: bytes) -> tuple[np.ndarray, int]:
    buf = bytearray(payload)
    root = encode.Get(packer.uoffset, buf, 0)
    matrix = Table(buf, root)
    n_rows = table_scalar(matrix, 4, number_types.Uint32Flags)
    n_cols = table_scalar(matrix, 6, number_types.Uint32Flags)
    columns_offset = number_types.UOffsetTFlags.py_type(matrix.Offset(8))
    if n_cols != 1 or matrix.VectorLen(columns_offset) != 1:
        raise RuntimeError(f"Expected one expression column; received {n_rows}x{n_cols}")

    column_pos = matrix.Indirect(matrix.Vector(columns_offset))
    column = Table(buf, column_pos)
    type_offset = number_types.UOffsetTFlags.py_type(column.Offset(4))
    union_type = column.Get(number_types.Uint8Flags, type_offset + column.Pos)
    if union_type != 1:  # NetEncoding.TypedArray.Float32Array
        raise RuntimeError(f"Expected Float32Array; received union type {union_type}")
    values = union_vector(column, 6, number_types.Float32Flags)

    index_type = table_scalar(matrix, 10, number_types.Uint8Flags)
    if index_type != 2:  # Int32Array
        raise RuntimeError(f"Expected Int32 column index; received union type {index_type}")
    col_index = union_vector(matrix, 12, number_types.Int32Flags)
    if len(values) != n_rows or len(col_index) != 1:
        raise RuntimeError("Malformed CELLxGENE matrix response")
    return values, int(col_index[0])


def fetch_gene(session: requests.Session, gene: str, force: bool = False) -> bytes:
    path = CACHE / f"{gene}.fbs"
    if not force and path.exists() and path.stat().st_size > 1_000_000:
        return path.read_bytes()
    response = session.get(
        URL,
        params={"var:var_names": gene},
        headers={"Accept": "application/octet-stream"},
        timeout=300,
    )
    response.raise_for_status()
    if len(response.content) < 1_000_000:
        raise RuntimeError(f"Unexpectedly short response for {gene}: {len(response.content)} bytes")
    path.write_bytes(response.content)
    return response.content


def main() -> None:
    metadata = pd.read_csv(META)
    rows = metadata["obs_index"].to_numpy(dtype=np.int64)
    CACHE.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0",
        # requests reliably decodes gzip/deflate.  Do not advertise Brotli
        # unless its optional decoder is installed.
        "Accept-Encoding": "gzip, deflate",
    })

    expression = metadata.copy()
    gene_indices: dict[str, int] = {}
    for i, gene in enumerate(CORE, start=1):
        payload = fetch_gene(session, gene)
        try:
            values, gene_index = decode_single_float_column(payload)
        except (RuntimeError, ValueError, IndexError, TypeError, ArithmeticError, struct.error):
            # Replace an incomplete or differently encoded cache atomically
            # through the normal downloader and validate the fresh response.
            payload = fetch_gene(session, gene, force=True)
            values, gene_index = decode_single_float_column(payload)
        if len(values) != N_ATLAS_CELLS:
            raise RuntimeError(f"Atlas row count changed for {gene}: {len(values)}")
        expression[gene] = np.asarray(values[rows], dtype=np.float32)
        gene_indices[gene] = gene_index
        print(f"[{i:02d}/{len(CORE)}] {gene}: {len(payload) / 2**20:.1f} MiB, var index {gene_index}")

    expression.to_csv(OUT, index=False, compression="gzip")
    print(f"Saved {len(expression):,} Schwann cells x {len(CORE)} genes to {OUT}")


if __name__ == "__main__":
    main()
