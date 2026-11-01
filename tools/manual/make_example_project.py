"""The example project of the manual: the synthetic export (make_example_data.py) set up as on the import page, and
analysed.

    python tools/manual/make_example_project.py [--no-run]

Folder: data/manual (gitignored). The project goes to data/manual/project/. Cellpose runs with the Python named in
data/manual/local.toml (`cellpose_python = '...'`, not in git: a path of this computer) or the environment variable
ISA_CELLPOSE_PYTHON; without either, with this environment's Python.
"""

from __future__ import annotations

import os
import sys
import tomllib
from pathlib import Path

from imagestream_analysis.config import CellOfInterest, Pool
from imagestream_analysis.project import DatasetEntry, Project, ZipEntry, run_project

DATA = Path(__file__).parents[2] / "data" / "manual"
FOLDER = DATA / "project"
CHANNELS = [1, 2, 3, 4, 7, 8]
NAMES = ["BF", "Marker", "Tumor", "Immune", "Subtype A", "Subtype B"]


def cellpose_python() -> str | None:
    local = DATA / "local.toml"
    if local.exists():
        with local.open("rb") as f:
            return str(tomllib.load(f).get("cellpose_python") or "") or None
    return os.environ.get("ISA_CELLPOSE_PYTHON")


def make() -> Project:
    p = Project(FOLDER)
    ch = p.settings.channels
    ch.brightfield, ch.marker_of_interest, ch.tumor = "BF", "Marker", "Tumor"
    ch.immune, ch.functional = ["Immune"], ["Subtype A", "Subtype B"]
    ch.cell_of_interest = CellOfInterest.TUMOR_CELL
    seg = p.settings.segmentation
    seg.diameter = 20.0
    python = cellpose_python()
    if python:
        seg.cellpose_python = Path(python)
    p.settings.classification.pool = Pool.EXPERIMENT
    p.settings.membrane.analysis = True
    p.settings.membrane.thickness = 3
    p.zips.append(ZipEntry(DATA / "example_export.zip", 0.5, list(CHANNELS), list(NAMES), {},
                           [DatasetEntry("Donor 1", "Example"), DatasetEntry("Donor 2", "Example")]))
    p.save()
    return p


if __name__ == "__main__":
    project = make()
    print(f"{FOLDER / 'project.toml'}; problems: {project.problems() or 'none'}")
    if "--no-run" not in sys.argv:
        print(run_project(FOLDER).to_string())
