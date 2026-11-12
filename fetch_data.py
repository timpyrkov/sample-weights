"""Download the NHANES table that correlations.py reads into data/ (not kept in git).

    uv run python fetch_data.py            # download and verify
    uv run python fetch_data.py --check    # verify only

Source: nhanes_userdata.csv.gz, one row per participant with named columns, sample weights,
strata and PSUs, published by pynhanes (https://pypi.org/project/pynhanes/) as its GitHub
release `data-v1` (built with pynhanes 1.1.0). The sha256 is fixed so every tester analyses
the same rows: a mismatch means the release changed. Report it, do not continue.
"""
from __future__ import annotations

import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"
TABLE = "nhanes_userdata.csv.gz"
URL = "https://github.com/timpyrkov/pynhanes/releases/download/data-v1/" + TABLE
SHA256 = "89e727a83bfeea89f20f8287027e2d19b8614de1ceca07dad18bf713bf4e4a11"


def problems() -> list[str]:
    """Empty if data/ holds the expected table, otherwise what is wrong."""
    path = DATA_DIR / TABLE
    if not path.exists():
        return [f"{TABLE}: missing"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256:
        return [f"{TABLE}: sha256 differs from the expected one"]
    return []


def main() -> None:
    if "--check" not in sys.argv:
        DATA_DIR.mkdir(exist_ok=True)
        print(f"downloading {TABLE} (about 30 MB)")
        tmp = DATA_DIR / (TABLE + ".part")
        with urllib.request.urlopen(URL, timeout=120) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f)
        tmp.replace(DATA_DIR / TABLE)
    bad = problems()
    print("\n".join(f"PROBLEM {b}" for b in bad) or "data/ is ready (sha256 verified)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
