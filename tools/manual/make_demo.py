"""The live demo of the manual: the dashboard of one dataset of the example project, as `manual-demo.zip` for the
`example-data` release (the site workflow unpacks it into manual/demo/), and in manual/demo/ for a local build.

    python tools/manual/make_demo.py

Run make_example_project.py first. The demo is not in git (manual/demo is ignored); upload it with
`gh release upload example-data data/manual/manual-demo.zip --clobber` BEFORE pushing, so that the site build of the
push takes the new demo.
"""

from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "manual" / "project" / "Example" / "Donor 1"
ZIP = ROOT / "data" / "manual" / "manual-demo.zip"
LOCAL = ROOT / "manual" / "demo"


def main() -> None:
    files = [SOURCE / "dashboard.html", *sorted((SOURCE / "dashboard_files").iterdir())]
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in files:
            z.write(f, f.relative_to(SOURCE).as_posix())
    if LOCAL.exists():
        shutil.rmtree(LOCAL)
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(LOCAL)
    print(f"{ZIP} ({ZIP.stat().st_size // 1024} kB, {len(files)} files); unpacked into {LOCAL}")


if __name__ == "__main__":
    main()
