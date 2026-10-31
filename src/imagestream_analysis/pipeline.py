"""The analysis of one sample, from the input stack to the output folder, in three stages:

1. measure: segmentation, clean-up and the per-cell measurements (`Measured`; the slow part);
2. classify: on its own (`run`), or together with the other datasets of a pool (project.run_project);
3. finish: interfaces, membranes and the output files (`Result`)."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import __version__, ideas
from .classification import class_map, classification_output, classify
from .config import Settings
from .dashboard import write_dashboard
from .interactions import add_neighbours, contact_table, interactions_per_event, multiplets, touching_pairs
from .interfaces import InterfaceMasks, doublets, event_table, interface_masks
from .io import EventStack, read_stack, save_stack, save_table
from .measure import edge_labels, measure_cells
from .membrane import (
    INTERFACE_OFFSET,
    MembraneResult,
    Profiles,
    analyze_membranes,
    expression_check,
    membrane_profiles,
    pair_of,
)
from .segmentation import clean_labels, markers, remove_labels, segment

log = logging.getLogger(__name__)

# table names of the macro (CELL_INTENSITY_TABLE, CELL_COUNT_TABLE, DOUBLETS_COUNT_TABLE)
CELL_TABLE = "Cell classification"
EVENT_TABLE = "All cell combinations"
DOUBLETS_TABLE = "Single Tumor Single Immune cell"  # macro v1.8: ... Single T cell
MEMBRANE_TABLE = "Membrane intensity analysis"
THRESHOLDS_TABLE = "Marker thresholds"
CONTACTS_TABLE = "Interactions"
MULTIPLETS_TABLE = "Multiplets"
AVERAGE_PROFILE = f"Average normalized membrane profile values - interfaces centered at {INTERFACE_OFFSET}"


@dataclass
class Measured:
    """The first stage: the cleaned labels and the measured cells, before any classification."""

    stack: EventStack
    labels: np.ndarray  # (T, Y, X) cleaned labels 1..N
    cells: pd.DataFrame  # measure.measure_cells
    cellpose: dict[str, Any]  # Cellpose version, device and timing ({} when the labels were loaded)


@dataclass
class Result:
    stack: EventStack
    labels: np.ndarray  # (T, Y, X) cleaned labels
    classes: np.ndarray  # (T, Y, X) 1 = immune cell, 2 = tumor cell, 3 = artefact
    cells: pd.DataFrame
    thresholds: pd.DataFrame  # the marker thresholds (of the pool)
    pairs: np.ndarray  # (P, 2) labels of the cells that touch
    contacts: pd.DataFrame  # interactions.contact_table
    multiplets: pd.DataFrame
    events: pd.DataFrame
    doublets: pd.DataFrame
    masks: InterfaceMasks
    out_dir: Path
    membranes: MembraneResult | None = None
    profiles: Profiles | None = None  # the membrane of every cell, for any pair of cell types (dashboard)
    expression: dict[str, Any] = field(default_factory=dict)  # where the marker of interest is (membrane.py)


def measure(stack: EventStack, settings: Settings, out_dir: Path) -> Measured:
    m = markers(stack, settings)
    raw, montage, info = segment(stack, m, settings)
    if montage is not None:  # as the macro: the Cellpose output, to reload later (segmentation.labelmap_file)
        save_stack(montage[None], out_dir / "labelmap_cells.tif", stack.pixel_size)
    seg = settings.segmentation
    labels = clean_labels(raw, seg.erosion, seg.min_cell_size, seg.max_cell_size)
    extents = stack.image_extents()
    edge = edge_labels(labels, extents, seg.erosion)
    if seg.exclude_edge_cells and edge.any():
        log.info("%d of %d cells are cut by the edge of their event image: left out", edge.sum(), len(edge))
        labels = remove_labels(labels, edge)
        edge = edge_labels(labels, extents, seg.erosion)
    else:
        log.info("%d of %d cells are cut by the edge of their event image (column Edge)", edge.sum(), len(edge))
    return Measured(stack, labels, measure_cells(labels, stack, m.tumor, m.immune, edge), info)


def analyze(stack: EventStack, settings: Settings, out_dir: Path) -> Result:
    """All stages for one dataset on its own."""
    measured = measure(stack, settings, out_dir)
    cells, thresholds = classify(measured.cells, settings.channels, settings.classification)
    return finish(measured, cells, thresholds, settings, out_dir)


def finish(measured: Measured, cells: pd.DataFrame, thresholds: pd.DataFrame, settings: Settings,
           out_dir: Path) -> Result:
    """The interfaces and membranes of the classified cells."""
    stack, labels = measured.stack, measured.labels
    pairs = touching_pairs(labels, settings.membrane.thickness)
    cells = add_neighbours(cells, pairs)
    classes = class_map(labels, cells)
    masks = interface_masks(classes, settings.membrane.thickness)
    moi = stack.channel(settings.channels.marker_of_interest)
    events = event_table(cells, stack.n_events, masks, moi)
    events["Artefacts"] = np.bincount(cells.loc[cells["Cell type"] == "Artefact", "frame"].to_numpy(dtype=int),
                                      minlength=stack.n_events)
    events["Interactions"] = interactions_per_event(cells, pairs, stack.n_events)
    result = Result(stack, labels, classes, cells, thresholds, pairs, contact_table(cells, pairs),
                    multiplets(cells, pairs, stack.n_events), events, doublets(events), masks, out_dir)
    ch = settings.channels
    if settings.membrane.analysis:
        m = settings.membrane
        result.profiles = membrane_profiles(labels, cells, moi, stack.median(ch.marker_of_interest),
                                            m.straightening_width, m.thickness, m.inside_only)
        result.membranes = analyze_membranes(cells, result.profiles, ch.cell_of_interest, m)
    wanted, partners = pair_of(settings.membrane, ch.cell_of_interest)
    moi_column = next(c for c in cells.columns if c.startswith(f"Mean ch{stack.number(ch.marker_of_interest)} ("))
    result.expression = expression_check(cells, result.profiles, moi_column, wanted, partners)
    level = logging.WARNING if result.expression["check"].startswith(("WARNING", "CAUTION")) else logging.INFO
    log.log(level, "Marker of interest, %s against %s: medians %.4g and %.4g (%s); %s", wanted, partners,
            result.expression["median (measured cells)"], result.expression["median (partners)"],
            result.expression["measured on"], result.expression["check"])
    result.cells.attrs["cellpose"] = measured.cellpose
    return result


def _with_object_ids(table: pd.DataFrame, ids: np.ndarray | None) -> pd.DataFrame:
    """Add the IDEAS object number of each row's event (column 'Frame', 1-based) at the end, when known."""
    if ids is None or "Frame" not in table:
        return table
    return table.assign(**{"Object ID": ids[table["Frame"].to_numpy(dtype=int) - 1]})


def save(result: Result, settings: Settings) -> None:
    out = result.out_dir
    ids = result.stack.object_ids
    save_table(_with_object_ids(classification_output(result.cells), ids), out / f"{CELL_TABLE}.tsv")
    save_table(result.thresholds, out / f"{THRESHOLDS_TABLE}.tsv")
    save_table(result.contacts, out / f"{CONTACTS_TABLE}.tsv")
    save_table(result.multiplets, out / f"{MULTIPLETS_TABLE}.tsv")
    save_table(_with_object_ids(result.events, ids), out / f"{EVENT_TABLE}.tsv")
    d = _with_object_ids(result.doublets, ids)
    d.insert(0, " ", np.arange(len(d)))  # the macro shows the row numbers in this table
    save_table(d, out / f"{DOUBLETS_TABLE}.tsv")
    save_stack(result.labels, out / "labelmap_final.tif", result.stack.pixel_size)
    save_stack(result.classes, out / "labelmap_classes.tif", result.stack.pixel_size)
    mem = result.membranes
    if mem is not None:
        save_table(_with_object_ids(mem.table, ids), out / f"{MEMBRANE_TABLE}.tsv")
        save_table(mem.average, out / f"{AVERAGE_PROFILE}.tsv")
        save_table(mem.histogram_interface, out / "Histogram_interface_stdDevs_above_mean.tsv")
        save_table(mem.histogram_control, out / "Histogram_control_stdDevs_above_mean.tsv")
        save_table(pd.DataFrame({"control stdDevs above mean": mem.control_scores}),
                   out / "Control windows stdDevs above mean.tsv")
        profiles = mem.profiles.astype(np.float32)[None]  # rows = the cells of the membrane table
        save_stack(profiles, out / "membrane_profiles_normalized.tif", result.stack.pixel_size)
    write_dashboard(result, settings)


def role_channels(settings: Settings) -> list[int]:
    """The channels with a role, by number (the settings must be resolved: ChannelSettings.resolved)."""
    return sorted({int(c) for refs in settings.channels.roles().values() for c in refs})


def load(settings: Settings) -> EventStack:
    """The events of the input: an IDEAS zip (the role channels plus input.channels) or a Combiner stack.
    Channels given by name in settings.channels are replaced by their numbers."""
    inp = settings.input
    path = Path(inp.image)
    if path.suffix.lower() != ".zip":
        stack = read_stack(path, inp.max_events)
        if inp.pixel_size:
            stack.pixel_size = inp.pixel_size
        settings.channels = settings.channels.resolved(stack.numbers, stack.channel_names)
        return stack
    if not inp.pixel_size:
        raise ValueError("input.pixel_size is required for an IDEAS zip (the files do not hold it)")
    # channels given by name: their IDEAS numbers, so that they are loaded
    settings.channels = settings.channels.resolved(inp.channels, inp.channel_names)
    datasets = ideas.scan(path)
    dataset = inp.dataset or (next(iter(datasets)) if len(datasets) == 1 else "")
    if not dataset:
        raise ValueError(f"{path.name} has {len(datasets)} datasets: set input.dataset ({', '.join(datasets)})")
    channels = sorted(set(inp.channels) | set(role_channels(settings)))
    name_of = dict(zip(inp.channels, inp.channel_names, strict=False))
    names = [name_of.get(c) or f"Ch{c}" for c in channels]
    return ideas.read_dataset(path, dataset, channels, names, inp.pixel_size, inp.max_events)


def start(settings: Settings) -> tuple[EventStack, Path]:
    """Load the input, make the output folder, start its log and save the settings."""
    stack = load(settings)
    out = Path(settings.input.output) / stack.name
    out.mkdir(parents=True, exist_ok=True)
    _setup_logging(out / "log.txt")
    log.info("imagestream-analysis %s: %s, %d of %d events, %d channels %s, pixel size %.4g µm", __version__,
             stack.name, stack.n_events, stack.n_events_in_file, len(stack.channel_names), stack.channel_names,
             stack.pixel_size)
    (out / "settings.toml").write_text(settings.to_toml(), encoding="utf-8")
    return stack, out


def run(settings: Settings) -> Result:
    """One dataset, classified on its own."""
    stack, out = start(settings)
    t0 = time.perf_counter()
    result = analyze(stack, settings, out)
    save(result, settings)
    write_run_info(result, settings, time.perf_counter() - t0)
    return result


def write_run_info(result: Result, settings: Settings, seconds: float, pooled_with: int = 0) -> None:
    cells = result.cells
    types = cells["Cell type"].value_counts()
    run_info = {"version": __version__, "input": str(settings.input.image), "events": result.stack.n_events,
                "cells": len(cells), "immune cells": int(types.get("Immune", 0)),
                "tumor cells": int(types.get("Tumor", 0)), "artefacts": int(types.get("Artefact", 0)),
                "subtypes": {k: int(v) for k, v in cells.loc[cells["Subtype"] != "", "Subtype"].value_counts().items()},
                "classified with": f"{pooled_with} other dataset(s)" if pooled_with else "this dataset only",
                "doublets": len(result.doublets), "marker of interest": result.expression,
                "membranes": None if result.membranes is None else {
                    "cells": len(result.membranes.table), "with_interface": int(result.membranes.has_interface.sum()),
                    "control_window_px": result.membranes.window},
                "cellpose": cells.attrs.get("cellpose", {}), "seconds": round(seconds, 1)}
    (result.out_dir / "run.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    log.info("Done in %.1f s: %d cells, %d single tumor + single immune cell events", run_info["seconds"],
             len(cells), len(result.doublets))


def _setup_logging(path: Path, append: bool = False) -> None:
    root = logging.getLogger()
    for h in list(root.handlers):
        if isinstance(h, logging.FileHandler):
            root.removeHandler(h)
            h.close()
    handler = logging.FileHandler(path, mode="a" if append else "w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    root.addHandler(handler)
    root.setLevel(logging.INFO)
