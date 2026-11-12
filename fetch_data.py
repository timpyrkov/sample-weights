"""Download the NHANES tables that correlations.py reads into data/ (not kept in git).

    uv run python fetch_data.py            # download and verify
    uv run python fetch_data.py --check    # verify only

Source: the parsed tables that pynhanes (https://pypi.org/project/pynhanes/) publishes as
its GitHub release `data-v1`, built with pynhanes 1.1.0:

- nhanes_userdata.csv.gz      one row per participant, named columns, weights, strata, PSUs
- nhanes_weights_dict.csv.gz  which sample weight belongs to each variable in each survey

Every tester must analyse the same rows, so the sha256 of each file is fixed below.
A mismatch means the release was rebuilt: report it, do not continue.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import urllib.request
from importlib import metadata
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"
RELEASE_URL = "https://github.com/timpyrkov/pynhanes/releases/download/data-v1/"
FILES = {  # name -> sha256
    "nhanes_userdata.csv.gz": "89e727a83bfeea89f20f8287027e2d19b8614de1ceca07dad18bf713bf4e4a11",
    "nhanes_weights_dict.csv.gz": "19f9d1a32042d5959ed0196beb158460bc87de2d30b5717a7594870983777fb8",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def problems() -> list[str]:
    """One message per file that is missing or differs from the expected sha256."""
    out = []
    for name, expected in FILES.items():
        path = DATA_DIR / name
        if not path.exists():
            out.append(f"{name}: missing")
        elif sha256(path) != expected:
            out.append(f"{name}: sha256 differs from the expected one")
    return out


def fetch() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    for name in FILES:
        print(f"downloading {name}")
        tmp = DATA_DIR / (name + ".part")
        with urllib.request.urlopen(RELEASE_URL + name, timeout=120) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f)
        tmp.replace(DATA_DIR / name)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="verify data/ without downloading")
    args = ap.parse_args()
    if not args.check:
        fetch()
    bad = problems()
    for line in bad:
        print("PROBLEM", line)
    if not bad:
        print(f"data/ is ready ({len(FILES)} files, sha256 verified, "
              f"installed pynhanes {metadata.version('pynhanes')})")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
