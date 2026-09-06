#!/usr/bin/env python3
"""Print targeted HDF5 metadata needed for CRC-Atlas expression extraction."""

from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location

import h5py


ROOT = Path(__file__).resolve().parents[1]
spec = spec_from_file_location("crc_remote", ROOT / "scripts/07_crc_atlas_remote_inspect.py")
module = module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

remote = module.HTTPRangeReader(
    "https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad",
    block_size=4 * 1024 * 1024,
    max_blocks=8,
)
with h5py.File(remote, "r") as h5:
    for path in ["X", "layers", "var", "var/var_names"]:
        obj = h5[path]
        print(path, type(obj).__name__, dict(obj.attrs))
        if isinstance(obj, h5py.Group):
            print(" keys", list(obj.keys()))
        else:
            print(" shape", obj.shape, "dtype", obj.dtype, "chunks", obj.chunks, "compression", obj.compression)
    if isinstance(h5["X"], h5py.Group):
        for key in h5["X"].keys():
            obj = h5["X"][key]
            print(f"X/{key}", obj.shape, obj.dtype, obj.chunks, obj.compression)
