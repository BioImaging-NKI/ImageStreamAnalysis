"""A project: several IDEAS zips (of one or more experiments) analysed with the same settings, and an overview of
their results. Stored as `project.toml` in the project's output folder; the results of a dataset go to
`<project folder>/<experiment>/<dataset>/`.

    [channels]                                           the panel: the role channels by name
    brightfield = 'BF'
    marker_of_interest = 'marker'
    tumor = 'tumor'
    immune = ['immune A', 'immune B']
    functional = ['subtype X']
    [segmentation] / [membrane] / [input]                the shared analysis settings (as a settings file)

    [[zip]]                                              one per IDEAS export
    path = 'D:/data/export.zip'
    pixel_size = 0.5
    channels = [1, 3, 6, 9, 11, 12]                       the IDEAS channels to load
    channel_names = ['BF', 'marker', 'immune A', 'tumor', 'immune B', 'subtype X']
    [[zip.dataset]]
    name = 'folder in the zip'
    experiment = 'label for the overview'
    include = true

The roles go by name, so a zip whose channels are numbered differently needs only its own names. A zip may still
have `roles = { tumor = 9, ... }` (channel numbers) that override the panel: projects written before 2026-10-31
gave every zip its roles by number; they are turned into a panel by name when all zips agree.
"""

from __future__ import annotations

import copy
import json
import logging
import re
import time
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import CellOfInterest, ChannelSettings, Settings, _from_dict, _toml_value
from .io import save_table

log = logging.getLogger(__name__)
PROJECT_FILE = "project.toml"
OVERVIEW_FILE = "overview.tsv"


@dataclass
class DatasetEntry:
    name: str
    experiment: str = ""
    include: bool = True


@dataclass
class ZipEntry:
    path: Path
    pixel_size: float | None = None
    channels: list[int] = field(default_factory=list)
    channel_names: list[str] = field(default_factory=list)
    roles: dict[str, Any] = field(default_factory=dict)  # overrides of the panel, by channel number (old projects)
    datasets: list[DatasetEntry] = field(default_factory=list)


@dataclass
class Project:
    folder: Path
    settings: Settings = field(default_factory=Settings)
    zips: list[ZipEntry] = field(default_factory=list)

    def samples(self) -> list[tuple[ZipEntry, DatasetEntry, Settings]]:
        """The included datasets with their complete settings."""
        out = []
        for z in self.zips:
            for d in z.datasets:
                if not d.include:
                    continue
                s = copy.deepcopy(self.settings)
                s.input.image = Path(z.path)
                s.input.dataset = d.name
                s.input.pixel_size = z.pixel_size
                s.input.channels = list(z.channels)
                s.input.channel_names = list(z.channel_names)
                s.input.output = self.folder / safe_name(d.experiment or Path(z.path).stem)
                s.channels = self.channels_of(z)
                out.append((z, d, s))
        return out

    def channels_of(self, z: ZipEntry) -> ChannelSettings:
        """The channel settings of a zip, by IDEAS channel number: the panel with the zip's overrides, the names
        looked up in the zip's channel names (ValueError if one is missing)."""
        ch = copy.deepcopy(self.settings.channels)
        for role, value in z.roles.items():
            setattr(ch, role, [int(v) for v in value] if isinstance(value, list) else int(value))
        return ch.resolved(z.channels, z.channel_names)

    def problems(self) -> list[str]:
        """What prevents a run (empty: ready)."""
        msgs = []
        if not any(d.include for z in self.zips for d in z.datasets):
            msgs.append("No datasets selected")
        if not self.settings.channels.immune:
            msgs.append("No immune-cell marker")
        for z in self.zips:
            name = Path(z.path).name
            if not z.pixel_size:
                msgs.append(f"{name}: pixel size not set")
            try:
                self.channels_of(z)
            except ValueError as e:
                msgs.append(f"{name}: {e}")
            if not Path(z.path).exists():
                msgs.append(f"{name}: file not found")
        return msgs

    def to_toml(self) -> str:
        lines = ["# ImageStream analysis project (imagestream-cli run-project)", "", self.settings.to_toml()]
        for z in self.zips:
            lines += ["[[zip]]", f"path = {_toml_value(str(z.path))}"]
            if z.pixel_size:
                lines.append(f"pixel_size = {z.pixel_size!r}")
            lines += [f"channels = {_toml_value(z.channels)}", f"channel_names = {_toml_value(z.channel_names)}"]
            if z.roles:
                lines.append("roles = { " + ", ".join(f"{k} = {_toml_value(v)}" for k, v in z.roles.items()) + " }")
            lines.append("")
            for d in z.datasets:
                lines += ["[[zip.dataset]]", f"name = {_toml_value(d.name)}",
                          f"experiment = {_toml_value(d.experiment)}", f"include = {_toml_value(d.include)}", ""]
        return "\n".join(lines)

    def save(self) -> Path:
        self.folder.mkdir(parents=True, exist_ok=True)
        path = self.folder / PROJECT_FILE
        path.write_text(self.to_toml(), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: str | Path) -> Project:
        path = Path(path)
        if path.is_dir():
            path = path / PROJECT_FILE
        with path.open("rb") as f:
            data = tomllib.load(f)
        return cls.from_dict(data, path.parent)

    @classmethod
    def from_dict(cls, data: dict[str, Any], folder: Path) -> Project:
        zips = [ZipEntry(Path(z["path"]), z.get("pixel_size"), [int(c) for c in z.get("channels", [])],
                         [str(n) for n in z.get("channel_names", [])],
                         ChannelSettings.migrate(dict(z.get("roles", {}))),
                         [DatasetEntry(d["name"], d.get("experiment", ""), bool(d.get("include", True)))
                          for d in z.get("dataset", [])])
                 for z in data.get("zip", [])]
        sections = {k: v for k, v in data.items() if k != "zip"}
        settings: Settings = _from_dict(Settings, sections)
        project = cls(folder, settings, zips)
        project.roles_to_panel()
        return project

    def roles_to_panel(self) -> None:
        """Turn per-zip roles by number (projects before 2026-10-31) into the panel by name, when every zip has the
        same named channel for each role; otherwise they stay as overrides."""
        if not self.zips or not all(z.roles for z in self.zips):
            return
        panels = []
        for z in self.zips:
            name_of = dict(zip(z.channels, z.channel_names, strict=False))
            panel: dict[str, Any] = {}
            for role, value in z.roles.items():
                names = [name_of.get(int(v), "") for v in (value if isinstance(value, list) else [value])]
                if not all(names):
                    return  # an unnamed role channel: keep the numbers
                panel[role] = names if role in ("immune", "functional") else names[0]
            panels.append(panel)
        if any(p != panels[0] for p in panels):
            return
        for role, value in panels[0].items():
            setattr(self.settings.channels, role, value)
        for z in self.zips:
            z.roles = {}


def safe_name(name: str) -> str:
    """A folder name: characters that Windows does not allow are replaced."""
    return re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", name).strip(" .") or "unnamed"


def summarize(result: Any, settings: Settings) -> dict[str, Any]:
    """Key numbers of one analysed dataset, for the overview."""
    cells = result.cells
    events = result.events
    d = result.doublets
    of_interest = "Tumor" if settings.channels.cell_of_interest == CellOfInterest.TUMOR_CELL else "Immune"
    side = "Tumor membrane ratio" if of_interest == "Tumor" else "Immune membrane ratio"
    row: dict[str, Any] = {
        "events": int(len(events)),
        "cells": int(len(cells)),
        "immune cells": int((cells["Cell type"] == "Immune").sum()),
        "tumor cells": int((cells["Cell type"] == "Tumor").sum()),
        "single tumor + single immune cell": int(len(d)),
        "doublet fraction": round(len(d) / len(events), 4) if len(events) else None,
        f"median {side} (doublets)": _median(d[side]) if len(d) else None,
    }
    mem = result.membranes
    if mem is not None:
        t = mem.table
        row |= {
            "membranes": int(len(t)),
            "with interface": int(mem.has_interface.sum()),
            "median stdDevs above mean": _median(t["stdDevs above mean"]),
            "median control (windows)": _median(pd.Series(mem.control_scores)),
        }
    return row


def _median(values: pd.Series) -> float | None:
    v = values.to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return round(float(np.median(v)), 4) if len(v) else None


def run_project(path: str | Path) -> pd.DataFrame:
    """Analyse all included datasets; a failing dataset is logged and skipped. Writes the overview."""
    from .pipeline import run  # scikit-learn and friends: only when running

    project = Project.load(path)
    problems = project.problems()
    if problems:
        raise ValueError("The project is not ready:\n- " + "\n- ".join(problems))
    samples = project.samples()
    for i, (z, d, settings) in enumerate(samples, start=1):
        print(f"ISA_PROGRESS {i}/{len(samples)} {d.experiment or Path(z.path).stem} / {d.name}", flush=True)
        t0 = time.perf_counter()
        try:
            result = run(settings)
        except Exception as e:  # noqa: BLE001 - one bad dataset must not stop the others
            log.exception("%s / %s failed: %s", Path(z.path).name, d.name, e)
            _write_info(Path(settings.input.output) / d.name, {"error": str(e)}, z, d)
            continue
        info = summarize(result, settings) | {"seconds": round(time.perf_counter() - t0, 1)}
        _write_info(result.out_dir, info, z, d)
    print(f"ISA_PROGRESS {len(samples)}/{len(samples)} done", flush=True)
    return write_overview(project.folder)


def _write_info(folder: Path, summary: dict[str, Any], z: ZipEntry, d: DatasetEntry) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    info = {"experiment": d.experiment, "zip": str(z.path), "dataset": d.name, "summary": summary}
    (folder / "summary.json").write_text(json.dumps(info, indent=2), encoding="utf-8")


def overview(folder: str | Path) -> pd.DataFrame:
    """One row per analysed dataset of the project (from the summary.json files of the result folders)."""
    folder = Path(folder)
    rows = []
    for f in sorted(folder.glob("*/*/summary.json")):
        info = json.loads(f.read_text(encoding="utf-8"))
        rows.append({"experiment": info.get("experiment", ""), "dataset": info.get("dataset", ""),
                     "zip": Path(info.get("zip", "")).name, "folder": str(f.parent.relative_to(folder))}
                    | info.get("summary", {}))
    return pd.DataFrame(rows)


def write_overview(folder: str | Path) -> pd.DataFrame:
    table = overview(folder)
    save_table(table, Path(folder) / OVERVIEW_FILE)
    return table
