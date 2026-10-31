"""The analysis of one sample, from the input stack to the output folder, in two stages:

1. measure: segmentation, clean-up and the per-cell measurements (`Measured`; the slow part);
2. the rest: classification, interfaces, membranes and the output files (`Result`)."""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import __version__, ideas
from .classification import class_map, classification_output, classify
from .config import Settings
from .dashboard import write_dashboard
from .interfaces import InterfaceMasks, doublets, event_table, interface_masks
from .io import EventStack, read_stack, save_stack, save_table
from .measure import edge_labels, measure_cells
from .membrane import INTERFACE_OFFSET, MembraneResult, analyze_membranes
from .segmentation import clean_labels, markers, remove_labels, segment

log = logging.getLogger(__name__)

# table names of the macro (CELL_INTENSITY_TABLE, CELL_COUNT_TABLE, DOUBLETS_COUNT_TABLE)
CELL_TABLE = "Cell classification"
EVENT_TABLE = "All cell combinations"
DOUBLETS_TABLE = "Single Tumor Single T cell"
MEMBRANE_TABLE = "Membrane intensity analysis"
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
    classes: np.ndarray  # (T, Y, X) 1 = T cell, 2 = tumor cell
    cells: pd.DataFrame
    events: pd.DataFrame
    doublets: pd.DataFrame
    masks: InterfaceMasks
    out_dir: Path
    membranes: MembraneResult | None = None


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
    return Measured(stack, labels, measure_cells(labels, stack, m.tumor, m.t_cell, edge), info)


def analyze(stack: EventStack, settings: Settings, out_dir: Path) -> Result:
    measured = measure(stack, settings, out_dir)
    labels = measured.labels
    cells = classify(measured.cells)
    classes = class_map(labels, cells)
    masks = interface_masks(classes, settings.membrane.thickness)
    moi = stack.channel(settings.channels.marker_of_interest)
    events = event_table(cells, stack.n_events, masks, moi)
    result = Result(stack, labels, classes, cells, events, doublets(events), masks, out_dir)
    ch = settings.channels
    if settings.membrane.analysis:
        interface = masks.t_cell_interface if ch.cell_of_interest == "T cell" else masks.tumor_interface
        result.membranes = analyze_membranes(labels, cells, interface, moi, stack.median(ch.marker_of_interest),
                                             ch.cell_of_interest, settings.membrane)
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
    ch = settings.channels
    return [ch.brightfield, ch.marker_of_interest, ch.tumor, ch.t_cell]


def load(settings: Settings) -> EventStack:
    """The events of the input: an IDEAS zip (the role channels plus input.channels) or a Combiner stack."""
    inp = settings.input
    path = Path(inp.image)
    if path.suffix.lower() != ".zip":
        stack = read_stack(path, inp.max_events)
        if inp.pixel_size:
            stack.pixel_size = inp.pixel_size
        return stack
    if not inp.pixel_size:
        raise ValueError("input.pixel_size is required for an IDEAS zip (the files do not hold it)")
    datasets = ideas.scan(path)
    dataset = inp.dataset or (next(iter(datasets)) if len(datasets) == 1 else "")
    if not dataset:
        raise ValueError(f"{path.name} has {len(datasets)} datasets: set input.dataset ({', '.join(datasets)})")
    channels = sorted(set(inp.channels) | set(role_channels(settings)))
    name_of = dict(zip(inp.channels, inp.channel_names, strict=False))
    names = [name_of.get(c) or f"Ch{c}" for c in channels]
    return ideas.read_dataset(path, dataset, channels, names, inp.pixel_size, inp.max_events)


def run(settings: Settings) -> Result:
    stack = load(settings)
    out = Path(settings.input.output) / stack.name
    out.mkdir(parents=True, exist_ok=True)
    _setup_logging(out / "log.txt")
    log.info("imagestream-analysis %s: %s, %d of %d events, %d channels %s, pixel size %.4g µm", __version__,
             stack.name, stack.n_events, stack.n_events_in_file, len(stack.channel_names), stack.channel_names,
             stack.pixel_size)
    (out / "settings.toml").write_text(settings.to_toml(), encoding="utf-8")
    t0 = time.perf_counter()
    result = analyze(stack, settings, out)
    save(result, settings)
    run_info = {"version": __version__, "input": str(settings.input.image), "events": stack.n_events,
                "cells": len(result.cells), "doublets": len(result.doublets),
                "membranes": None if result.membranes is None else {
                    "cells": len(result.membranes.table), "with_interface": int(result.membranes.has_interface.sum()),
                    "control_window_px": result.membranes.window},
                "cellpose": result.cells.attrs.get("cellpose", {}), "seconds": round(time.perf_counter() - t0, 1)}
    (out / "run.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    log.info("Done in %.1f s: %d cells, %d single tumor + single T cell events", run_info["seconds"],
             len(result.cells), len(result.doublets))
    return result


def _setup_logging(path: Path) -> None:
    root = logging.getLogger()
    for h in list(root.handlers):
        if isinstance(h, logging.FileHandler):
            root.removeHandler(h)
            h.close()
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    root.addHandler(handler)
    root.setLevel(logging.INFO)
