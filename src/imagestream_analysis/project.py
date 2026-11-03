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
    condition = 'treated'                                optional: the comparison of conditions
    replicate = 'donor 1'                                (comparison.py); [comparison] reference = 'control'

The roles go by name, so a zip whose channels are numbered differently needs only its own names. A zip may still
have `roles = { tumor = 9, ... }` (channel numbers) that override the panel: projects written before 2026-10-31
gave every zip its roles by number; they are turned into a panel by name when all zips agree.
"""

from __future__ import annotations

import copy
import json
import logging
import math
import re
import time
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .comparison import write_comparison
from .config import CellOfInterest, ChannelSettings, Pool, Settings, _from_dict, _toml_value
from .io import save_table
from .membrane import pair_of

log = logging.getLogger(__name__)
PROJECT_FILE = "project.toml"
OVERVIEW_FILE = "overview.tsv"
INTERACTIONS_FILE = "interactions.tsv"
SUBSETS_FILE = "subset preference.tsv"
PAIRS_FILE = "membrane pairs.tsv"


@dataclass
class DatasetEntry:
    name: str
    experiment: str = ""
    include: bool = True
    condition: str = ""  # control, treated, ...: compared across the replicates (comparison.py)
    replicate: str = ""  # the independent unit: a donor, an independent experiment


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
    reference: str = ""  # the condition the others are compared with (empty: the first one)

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
        if self.settings.membrane.analysis and not self.settings.channels.has_marker_of_interest:
            msgs.append("The membrane analysis needs a marker of interest (or switch it off)")
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
        if self.reference:
            lines += ["[comparison]", f"reference = {_toml_value(self.reference)}", ""]
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
                          f"experiment = {_toml_value(d.experiment)}", f"include = {_toml_value(d.include)}"]
                lines += [f"{k} = {_toml_value(v)}" for k, v in (("condition", d.condition),
                                                                  ("replicate", d.replicate)) if v]
                lines.append("")
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
                         [DatasetEntry(d["name"], d.get("experiment", ""), bool(d.get("include", True)),
                                       str(d.get("condition", "")).strip(), str(d.get("replicate", "")).strip())
                          for d in z.get("dataset", [])])
                 for z in data.get("zip", [])]
        sections = {k: v for k, v in data.items() if k not in ("zip", "comparison")}
        settings: Settings = _from_dict(Settings, sections)
        project = cls(folder, settings, zips, str(data.get("comparison", {}).get("reference", "")).strip())
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
        "artefacts": int((cells["Cell type"] == "Artefact").sum()),
        "single tumor + single immune cell": int(len(d)),
        "interactions": int(len(result.pairs)),
        "MOI partners / measured": _finite(result.expression.get("partners / measured")),
        "marker check": result.expression.get("level", ""),
        "tumor–immune contacts": int(result.contacts.loc[(result.contacts["cell type"] == "Tumor")
                                                         & result.contacts["neighbour type"].str.startswith("Immune"),
                                                         "contacts"].sum()),
        "doublet fraction": round(len(d) / len(events), 4) if len(events) else None,
    }
    if side in d:
        row[f"median {side} (doublets)"] = _median(d[side]) if len(d) else None
    for subtype, count in cells.loc[cells["Subtype"] != "", "Subtype"].value_counts().items():
        row[f"immune {subtype}"] = int(count)
    mem = result.membranes
    if mem is not None:
        t = mem.table
        cells_of, partners = pair_of(settings.membrane, settings.channels.cell_of_interest)
        row |= {
            "membrane pair": f"{' / '.join(cells_of)} with {' / '.join(partners)}",
            "membranes": int(len(t)),
            "with interface": int(mem.has_interface.sum()),
            "median stdDevs above mean": _median(t["stdDevs above mean"]),
            "median control (windows)": _median(pd.Series(mem.control_scores)),
        }
    return row


def _finite(value: Any) -> float | None:
    """A number for the overview: rounded, or None when there is none (missing, NaN or infinite)."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return round(v, 3) if math.isfinite(v) else None


def _median(values: pd.Series) -> float | None:
    v = values.to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return round(float(np.median(v)), 4) if len(v) else None


def run_project(path: str | Path, only: list[str] | None = None) -> pd.DataFrame:
    """Analyse all included datasets; a failing dataset is logged and skipped. Writes the overview.

    classification.pool = "experiment": all datasets are measured first, then the cells of each experiment are
    classified together (one set of thresholds, one clustering), then every dataset is finished. "dataset": each
    dataset is analysed on its own.

    only: the names of the datasets to run (again); the others keep their results, and the overview is written from
    all of them. With pools per experiment, a dataset is classified with the other datasets of its experiment, so
    its whole experiment runs."""
    project = Project.load(path)
    problems = project.problems()
    if problems:
        raise ValueError("The project is not ready:\n- " + "\n- ".join(problems))
    samples = project.samples()
    if only:
        samples = select_samples(samples, only, project.settings.classification.pool)
    if project.settings.classification.pool == Pool.DATASET:
        _run_each(samples)
    else:
        _run_pooled(samples)
    print(f"ISA_PROGRESS {len(samples)}/{len(samples)} done", flush=True)
    return write_overview(project.folder)


Sample = tuple["ZipEntry", "DatasetEntry", Settings]


def select_samples(samples: list[Sample], only: list[str], pool: Pool) -> list[Sample]:
    """The samples of the named datasets; with pools per experiment, all datasets of their experiments."""
    names = set(only)
    chosen = [s for s in samples if s[1].name in names]
    missing = names - {d.name for _, d, _ in chosen}
    if missing:
        raise ValueError(f"No included dataset named {', '.join(sorted(missing))}")
    if pool == Pool.EXPERIMENT:
        experiments = {d.experiment or Path(z.path).stem for z, d, _ in chosen}
        chosen = [s for s in samples if (s[1].experiment or Path(s[0].path).stem) in experiments]
        extra = [d.name for _, d, _ in chosen if d.name not in names]
        if extra:
            log.info("Classified together with the other datasets of the experiment, so these run too: %s",
                     ", ".join(extra))
    return chosen


def _label(z: ZipEntry, d: DatasetEntry) -> str:
    return f"{d.experiment or Path(z.path).stem} / {d.name}"


def _failed(z: ZipEntry, d: DatasetEntry, settings: Settings, e: Exception) -> None:
    log.exception("%s / %s failed: %s", Path(z.path).name, d.name, e)
    _write_info(Path(settings.input.output) / d.name, {"error": str(e)}, z, d)


def _run_each(samples: list[Sample]) -> None:
    from .pipeline import run  # scikit-learn and friends: only when running

    for i, (z, d, settings) in enumerate(samples, start=1):
        print(f"ISA_PROGRESS {i - 1}/{len(samples)} {_label(z, d)}", flush=True)
        t0 = time.perf_counter()
        try:
            result = run(settings)
        except Exception as e:  # noqa: BLE001 - one bad dataset must not stop the others
            _failed(z, d, settings, e)
            continue
        _write_info(result.out_dir, summarize(result, settings) | {"seconds": round(time.perf_counter() - t0, 1)},
                    z, d)


def _run_pooled(samples: list[Sample]) -> None:
    from .classification import classify_pool
    from .pipeline import Measured, _setup_logging, finish, load, measure, save, start, write_run_info

    n = len(samples)
    measured: dict[int, tuple[Measured, Path, float]] = {}
    for i, (z, d, settings) in enumerate(samples):  # 1. segment and measure every dataset
        print(f"ISA_PROGRESS {i}/{n} measuring {_label(z, d)}", flush=True)
        t0 = time.perf_counter()
        try:
            stack, out = start(settings)
            m = measure(stack, settings, out)
        except Exception as e:  # noqa: BLE001 - one bad dataset must not stop the others
            _failed(z, d, settings, e)
            continue
        m.stack.data = m.stack.data[:0]  # the images are read again to finish: keep the memory free meanwhile
        measured[i] = (m, out, time.perf_counter() - t0)
    groups: dict[str, list[int]] = {}
    for i in measured:
        z, d, _ = samples[i]
        groups.setdefault(d.experiment or Path(z.path).stem, []).append(i)
    for experiment, members in groups.items():  # 2. classify each experiment together, 3. finish each dataset
        t0 = time.perf_counter()
        try:
            pooled, thresholds = classify_pool([(measured[i][0].cells, samples[i][2].channels) for i in members],
                                               samples[members[0]][2].classification)
        except Exception as e:  # noqa: BLE001 - a bad experiment must not stop the others
            for i in members:
                _failed(samples[i][0], samples[i][1], samples[i][2], e)
            continue
        log.info("Experiment '%s': %d dataset(s) classified together", experiment, len(members))
        share = (time.perf_counter() - t0) / len(members)
        for i, cells in zip(members, pooled, strict=True):
            z, d, settings = samples[i]
            m, out, seconds = measured[i]
            print(f"ISA_PROGRESS {i}/{n} finishing {_label(z, d)}", flush=True)
            t1 = time.perf_counter()
            try:
                _setup_logging(out / "log.txt", append=True)
                m.stack = load(copy.deepcopy(settings))  # the images again
                cells.attrs["cellpose"] = m.cellpose
                result = finish(m, cells, thresholds, settings, out)
                save(result, settings)
                total = seconds + share + time.perf_counter() - t1
                write_run_info(result, settings, total, pooled_with=len(members) - 1)
            except Exception as e:  # noqa: BLE001 - one bad dataset must not stop the others
                _failed(z, d, settings, e)
                continue
            _write_info(out, summarize(result, settings) | {"seconds": round(total, 1)}, z, d)


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
    table = pd.DataFrame(rows)
    # the immune subtypes: 0 where a dataset has none, in one order (negative, then by number of markers, by name)
    subtypes = [c for c in table.columns if c.startswith("immune ") and c != "immune cells"]
    if subtypes:
        table[subtypes] = table[subtypes].fillna(0).astype(int)
        ordered = sorted(subtypes, key=lambda c: (c != "immune negative", c.count("+"), c))
        rest = [c for c in table.columns if c not in subtypes]
        at = list(table.columns).index(subtypes[0])
        at = len([c for c in table.columns[:at] if c not in subtypes])
        table = table[rest[:at] + ordered + rest[at:]]
    return table


def write_overview(folder: str | Path) -> pd.DataFrame:
    """overview.tsv, and interactions.tsv: the contact tables of all datasets, and summed per experiment;
    subset preference.tsv, membrane pairs.tsv, and condition comparison.tsv (with a project file)."""
    folder = Path(folder)
    table = overview(folder)
    save_table(table, folder / OVERVIEW_FILE)
    parts = []
    for row in table.itertuples():
        f = folder / str(row.folder) / "Interactions.tsv"
        if f.exists():
            parts.append(pd.read_csv(f, sep="	").assign(experiment=row.experiment, dataset=row.dataset))
    if parts:
        each = pd.concat(parts, ignore_index=True)
        keys = ["experiment", "cell type", "neighbour type"]
        summed = (each.groupby(keys, as_index=False)[["contacts", "cells touching", "contacts (3+ cells)",
                                                      "expected (3+ cells)"]].sum()
                  .merge(each.drop_duplicates(["experiment", "dataset", "cell type"])
                         .groupby(["experiment", "cell type"], as_index=False)["cells"].sum(), on=keys[:2]))
        summed["fraction touching"] = summed["cells touching"] / summed["cells"]
        summed["observed / expected"] = (summed["contacts (3+ cells)"]
                                         / summed["expected (3+ cells)"].where(summed["expected (3+ cells)"] > 0))
        summed.insert(1, "dataset", "(all)")
        columns = ["experiment", "dataset", "cell type", "neighbour type", "cells", "contacts", "cells touching",
                   "fraction touching", "contacts (3+ cells)", "expected (3+ cells)", "observed / expected"]
        save_table(pd.concat([each[columns], summed[columns]], ignore_index=True), folder / INTERACTIONS_FILE)
    write_subsets(folder, table)
    pairs = [pd.read_csv(f, sep="\t").assign(experiment=row.experiment, dataset=row.dataset)
             for row in table.itertuples() if (f := folder / str(row.folder) / "Membrane pairs.tsv").exists()]
    if pairs:
        each = pd.concat(pairs, ignore_index=True)
        save_table(each[["experiment", "dataset", *[c for c in each.columns if c not in ("experiment", "dataset")]]],
                   folder / PAIRS_FILE)
    if (folder / PROJECT_FILE).exists() and len(table):
        write_comparison(Project.load(folder), table)
    return table


def write_subsets(folder: Path, table: pd.DataFrame) -> None:
    """subset preference.tsv: the subset preference of every dataset (with its 95 % interval), and per experiment the
    observed and expected contacts summed over its datasets (their ratio; no interval: compare the datasets)."""
    parts = []
    for row in table.itertuples():
        f = folder / str(row.folder) / "Subset preference.tsv"
        if f.exists():
            parts.append(pd.read_csv(f, sep="\t").assign(experiment=row.experiment, dataset=row.dataset))
    if not parts:
        return
    each = pd.concat(parts, ignore_index=True)
    summed = each.groupby(["experiment", "cell type", "neighbour type"], as_index=False)[
        ["contacts", "expected (subtypes shuffled)"]].sum()
    exp = summed["expected (subtypes shuffled)"]
    summed["observed / expected (subtypes shuffled)"] = summed["contacts"] / exp.where(exp > 0)
    summed.insert(1, "dataset", "(all)")
    columns = ["experiment", "dataset", "cell type", "neighbour type", "contacts", "expected (subtypes shuffled)",
               "observed / expected (subtypes shuffled)", "95% low", "95% high"]
    save_table(pd.concat([each, summed], ignore_index=True).reindex(columns=columns), folder / SUBSETS_FILE)
