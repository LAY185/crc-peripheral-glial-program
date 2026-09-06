from __future__ import annotations

import importlib
import json
import platform
import sys
from pathlib import Path


REQUIRED = [
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "seaborn",
    "statsmodels",
    "h5py",
    "flatbuffers",
    "requests",
    "yaml",
]


def main() -> None:
    versions: dict[str, str] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    missing: list[str] = []
    for name in REQUIRED:
        try:
            module = importlib.import_module(name)
            versions[name] = getattr(module, "__version__", "available")
        except Exception:
            missing.append(name)

    for directory in [
        "data/raw",
        "data/processed",
        "results/figures",
        "results/tables",
        "results/cache",
        "manuscript",
        "logs",
    ]:
        Path(directory).mkdir(parents=True, exist_ok=True)

    report = {"versions": versions, "missing": missing}
    Path("logs/environment.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
