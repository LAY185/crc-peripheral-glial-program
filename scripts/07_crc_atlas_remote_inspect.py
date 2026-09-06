#!/usr/bin/env python3
"""Inspect the public CRC Atlas h5ad through HTTP byte-range requests.

This avoids downloading the 32.7 GB atlas merely to determine whether its
Schwann-cell compartment is useful for the project.
"""

import io
from collections import OrderedDict
from urllib.request import Request, urlopen

import h5py


URL = "https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad"


class HTTPRangeReader(io.RawIOBase):
    """Small seekable HTTP file backed by aligned Range requests."""

    def __init__(self, url, block_size=8 * 1024 * 1024, max_blocks=8):
        self.url = url
        self.block_size = block_size
        self.max_blocks = max_blocks
        self.pos = 0
        req = Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=60) as response:
            self.size = int(response.headers["Content-Length"])
        self.cache = OrderedDict()

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        if whence == io.SEEK_SET:
            self.pos = offset
        elif whence == io.SEEK_CUR:
            self.pos += offset
        elif whence == io.SEEK_END:
            self.pos = self.size + offset
        else:
            raise ValueError("invalid whence")
        return self.pos

    def _block(self, index):
        if index in self.cache:
            self.cache.move_to_end(index)
            return self.cache[index]
        start = index * self.block_size
        end = min(start + self.block_size, self.size) - 1
        req = Request(
            self.url,
            headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"},
        )
        with urlopen(req, timeout=120) as response:
            payload = response.read()
        expected = end - start + 1
        if len(payload) != expected:
            raise IOError(f"range {start}-{end}: expected {expected}, got {len(payload)}")
        self.cache[index] = payload
        self.cache.move_to_end(index)
        while len(self.cache) > self.max_blocks:
            self.cache.popitem(last=False)
        return payload

    def read(self, size=-1):
        if size is None or size < 0:
            size = self.size - self.pos
        size = min(size, self.size - self.pos)
        if size <= 0:
            return b""
        chunks = []
        remaining = size
        while remaining:
            index = self.pos // self.block_size
            block = self._block(index)
            within = self.pos % self.block_size
            take = min(remaining, len(block) - within)
            chunks.append(block[within:within + take])
            self.pos += take
            remaining -= take
        return b"".join(chunks)

    def readinto(self, b):
        payload = self.read(len(b))
        b[:len(payload)] = payload
        return len(payload)


def describe(name, obj):
    if isinstance(obj, h5py.Dataset):
        print("DATASET", name, obj.shape, obj.dtype, dict(obj.attrs))
    else:
        print("GROUP", name, dict(obj.attrs))


def main():
    remote = HTTPRangeReader(URL)
    with h5py.File(remote, "r") as h5:
        print("ROOT", list(h5.keys()), dict(h5.attrs))
        for key in ["X", "obs", "var", "raw"]:
            if key in h5:
                print("\n##", key)
                h5[key].visititems(describe)


if __name__ == "__main__":
    main()
