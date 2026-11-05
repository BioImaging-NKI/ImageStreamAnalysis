"""The analysis of one sample, from the input stack to the output folder, in three stages:

1. measure: segmentation, clean-up and the per-cell measurements (`Measured`; the slow part);
2. classify: on its own (`run`), or together with the other datasets of a pool (project.run_project);
3. finish: interfaces, membranes (for every marker of interest) and the output files (`Result`)."""

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
from .config import Settings, membrane_px, straightening_px
from .dashboard import write_dashboard
from .interactions import (
    add_neighbours,
    contact_table,
    interactions_per_event,
    multiplets,
    subset_preference,
    touching_pairs,
)
from .interfaces import InterfaceMasks, doublets, event_table, interface_masks
from .io import EventStack, read_stack, safe_name, save_stack, save_table
from .measure import edge_labels, measure_cells
from .membrane import (
    INTERFACE_OFFSET,
    MembraneResult,
    Profiles,
    analyze_membranes,
    carrier,
    expression_check,
    location_notes,
    membrane_profiles_of,
    pair_of,
    pair_table,
    select,
    selector_problems,
)
from .segmentation import clean_labels, markers, remove_labels, save_cellpose_result, segment

log = logging.getLogger(__name__)

# table names of the macro (CELL_INTENSITY_TABLE, CELL_COUNT_TABLE, DOUBLETS_COUNT_TABLE)
CELL_TABLE = "Cell classification"
EVENT_TABLE = "All cell combinations"
DOUBLETS_TABLE = "Single Tumor Single Immune cell"  # macro v1.8: ... Single T cell
MEMBRANE_TABLE = "Membrane intensity analysis"
LOCATION_TABLE = "Membrane location"  # where the membrane lies across the outline (membrane.Location)
THRESHOLDS_TABLE = "Marker thresholds"
CONTACTS_TABLE = "Interactions"
SUBSETS_TABLE = "Subset preference"
PAIRS_TABLE = "Membrane pairs"
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
class MarkerResult:
    """One marker of interest: the cell type that carries it, the pair measured, the marker check and the membranes."""

    name: str  # the channel name
    number: int  # the channel number
    column: str  # its column in the cell table, "Mean ch<N> (<name>)"
    carried_by: str  # "Tumor" or "Immune" ("" when undecided): the type with the higher median (membrane.carrier)
    medians: dict[str, float]  # that median, of the tumor and of the immune cells
    cells: list[str]  # the pair measured (membrane.pair_of)
    partners: list[str]
    expression: dict[str, Any]  # the marker check of the pair (membrane.expression_check)
    profiles: Profiles | None = None  # the membrane of every cell, for any pair of cell types (dashboard)
    membranes: MembraneResult | None = None
    pairs_table: pd.DataFrame = field(default_factory=pd.DataFrame)  # membrane.pair_table
    suffix: str = ""  # of its output files: " - <name>" (every marker alike, decided 2026-11-03)


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
    subsets: pd.DataFrame = field(default_factory=pd.DataFrame)  # interactions.subset_preference
    markers: list[MarkerResult] = field(default_factory=list)  # one per marker of interest, in the settings' order

    @property
    def pairs_table(self) -> pd.DataFrame:
        """membrane.pair_table of every marker, with the columns marker and marker on in front (Membrane pairs.tsv)."""
        parts = [m.pairs_table.assign(**{"marker": m.name, "marker on": m.carried_by})
                 for m in self.markers if len(m.pairs_table)]
        if not parts:
            return pd.DataFrame()
        t = pd.concat(parts, ignore_index=True)
        return t[["marker", "marker on", *[c for c in t.columns if c not in ("marker", "marker on")]]]


def measure(stack: EventStack, settings: Settings, out_dir: Path) -> Measured:
    m = markers(stack, settings)
    raw, montage, info = segment(stack, m, settings, previous=out_dir)
    if montage is not None:  # the Cellpose output: reused by the next run, or loaded with segmentation.labelmap_file
        save_cellpose_result(montage, info, out_dir, stack.pixel_size)
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
    reach = membrane_px(settings.membrane.thickness, stack.pixel_size)
    pairs = touching_pairs(labels, reach)
    cells = add_neighbours(cells, pairs)
    classes = class_map(labels, cells)
    masks = interface_masks(classes, reach)
    refs = list(settings.channels.marker_of_interest)
    if not refs and settings.membrane.analysis:
        raise ValueError("the membrane analysis needs a marker of interest (channels.marker_of_interest)")
    mois = [stack.channel(r) for r in refs]
    named = marker_names(cells, [stack.number(r) for r in refs])
    events = event_table(cells, stack.n_events, masks, [(name, moi) for (_, _, name), moi in zip(named, mois,
                                                                                              strict=True)])
    events["Artefacts"] = np.bincount(cells.loc[cells["Cell type"] == "Artefact", "frame"].to_numpy(dtype=int),
                                      minlength=stack.n_events)
    events["Interactions"] = interactions_per_event(cells, pairs, stack.n_events)
    result = Result(stack, labels, classes, cells, thresholds, pairs, contact_table(cells, pairs),
                    multiplets(cells, pairs, stack.n_events), events, doublets(events), masks, out_dir,
                    subsets=subset_preference(cells, pairs))
    result.cells.attrs["cellpose"] = measured.cellpose
    if not refs:
        log.info("No marker of interest: cell types and contacts only")
        return result
    m = settings.membrane
    profiles: list[Profiles | None] = [None] * len(refs)
    if m.analysis:  # the outlines and neighbour tracks once, a profile per marker
        width = straightening_px(m.straightening_width, stack.pixel_size)
        profiles = list(membrane_profiles_of(labels, cells, mois, [stack.median(r) for r in refs], width, reach,
                                             m.inside_only, str(m.band_centre), [n for _, _, n in named],
                                             stack.pixel_size))
        location = profiles[0].location if profiles and profiles[0] is not None else None
        if location is not None:
            for level, text in location_notes(location, width, stack.pixel_size):
                log.log(level, text)
    for (number, column, name), prof in zip(named, profiles, strict=True):
        result.markers.append(marker_result(cells, number, column, name, prof, settings))
    return result


def marker_names(cells: pd.DataFrame, numbers: list[int]) -> list[tuple[int, str, str]]:
    """(channel number, cell-table column "Mean ch<N> (<name>)", name) of every marker of interest; a name that two
    markers share gets its channel number added, so that columns and files stay apart."""
    out = []
    for number in numbers:
        prefix = f"Mean ch{number} ("
        column = next(c for c in cells.columns if c.startswith(prefix))
        out.append((number, column, column[len(prefix) : -1]))
    names = [n for _, _, n in out]
    return [(k, c, n if names.count(n) == 1 else f"{n} (ch{k})") for k, c, n in out]


CARRIER_CLOSE = 0.8  # the other type at this fraction of the carrier's level or more: a warning (chosen narrowly)


def marker_result(cells: pd.DataFrame, number: int, column: str, name: str, profiles: Profiles | None,
                  settings: Settings) -> MarkerResult:
    """One marker of interest: the cell type that carries it (unless membrane.cells is set), its pair, the marker
    check and, with profiles, the membrane analysis."""
    suffix = f" - {safe_name(name)}"
    m = settings.membrane
    # the type that carries the marker: with membrane profiles, the type whose free membrane is the brightest in
    # the marker averaged across the outline (it does not depend on the band's width or place, 2026-11-05: with the
    # band's own free-membrane values it flipped with those settings); else the whole-cell medians
    location = profiles.location if profiles is not None else None
    carried_by, medians = location.carrier(name) if location is not None else carrier(cells, None, column)
    if not carried_by and location is not None:
        carried_by, medians = carrier(cells, None, column)
    wanted, partners = pair_of(m, carried_by, name)
    source = ("height of the free membrane across the outline" if location is not None
              and location.carrier(name)[0] else "median of the whole cell")
    values = sorted((v for v in medians.values() if np.isfinite(v)), reverse=True)
    if len(values) == 2 and values[0] > 0 and values[1] >= CARRIER_CLOSE * values[0]:
        log.warning("%s: the tumor and immune cells carry it at a similar level (%s; %s): the measured cells (%s) "
                    "are chosen by a small margin. Set the pair (membrane.pairs) to fix it.", name,
                    ", ".join(f"{k.lower()} {v:.4g}" for k, v in medians.items()), source, carried_by.lower())
    set_cells, set_partners = m.pair_settings(name)
    chosen = ("set in membrane.pairs" if name in m.pairs else "set in membrane.cells / partners"
              if set_cells or set_partners else "automatic: the type that carries it")
    log.info("%s: carried by the %s cells (%s; %s): membranes of %s with %s (%s)", name,
             carried_by.lower() or "(no tumor or immune)",
             ", ".join(f"{k.lower()} {v:.4g}" for k, v in medians.items()), source, " / ".join(wanted),
             " / ".join(partners), chosen)
    expression = expression_check(cells, profiles, column, wanted, partners)
    level = (logging.WARNING if expression["check"].startswith(("WARNING", "CAUTION", "NOT INFORMATIVE"))
             else logging.INFO)
    log.log(level, "%s, %s against %s: medians %.4g and %.4g (%s); %s", name, wanted, partners,
            expression["median (measured cells)"], expression["median (partners)"], expression["measured on"],
            expression["check"])
    out = MarkerResult(name, number, column, carried_by, medians, wanted, partners,
                       {"marker": name, "marker on": carried_by} | expression, profiles, suffix=suffix)
    for what, sels in (("cells", wanted), ("partners", partners)):  # a valid subset may have no cells here
        if not select(cells, sels).any():
            log.warning("%s: no %s of the pair in this dataset (%s): no membrane results for this marker", name,
                        what, " / ".join(sels))
    if profiles is not None:
        out.membranes = analyze_membranes(cells, profiles, carried_by, m, marker=name)
        out.pairs_table = pair_table(cells, profiles, column, (wanted, partners))
    return out


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
    if len(result.subsets):
        save_table(result.subsets, out / f"{SUBSETS_TABLE}.tsv")
    if len(result.pairs_table):
        save_table(result.pairs_table, out / f"{PAIRS_TABLE}.tsv")
    save_table(result.multiplets, out / f"{MULTIPLETS_TABLE}.tsv")
    save_table(_with_object_ids(result.events, ids), out / f"{EVENT_TABLE}.tsv")
    d = _with_object_ids(result.doublets, ids)
    d.insert(0, " ", np.arange(len(d)))  # the macro shows the row numbers in this table
    save_table(d, out / f"{DOUBLETS_TABLE}.tsv")
    save_stack(result.labels, out / "labelmap_final.tif", result.stack.pixel_size)
    save_stack(result.classes, out / "labelmap_classes.tif", result.stack.pixel_size)
    for marker in result.markers:  # the macro's file names, with " - <marker>" added
        mem, sfx = marker.membranes, marker.suffix
        if mem is None:
            continue
        if not len(mem.table):  # no measured cells: no (empty) files, and none left from an earlier run
            log.warning("%s: no measured cells (membrane of %s): no membrane files", marker.name,
                        " / ".join(marker.cells))
            for f in out.iterdir():
                if f.stem.endswith(sfx) and f.name.startswith((MEMBRANE_TABLE, AVERAGE_PROFILE, "Histogram_",
                                                               "Control windows", "membrane_profiles_normalized")):
                    f.unlink()
            continue
        save_table(_with_object_ids(mem.table, ids), out / f"{MEMBRANE_TABLE}{sfx}.tsv")
        save_table(mem.average, out / f"{AVERAGE_PROFILE}{sfx}.tsv")
        save_table(mem.histogram_interface, out / f"Histogram_interface_stdDevs_above_mean{sfx}.tsv")
        save_table(mem.histogram_control, out / f"Histogram_control_stdDevs_above_mean{sfx}.tsv")
        save_table(pd.DataFrame({"control stdDevs above mean": mem.control_scores}),
                   out / f"Control windows stdDevs above mean{sfx}.tsv")
        profiles = mem.profiles.astype(np.float32)[None]  # rows = the cells of the membrane table
        save_stack(profiles, out / f"membrane_profiles_normalized{sfx}.tif", result.stack.pixel_size)
    location = next((mk.profiles.location for mk in result.markers if mk.profiles is not None
                     and mk.profiles.location is not None), None)
    if location is not None:  # where the membrane lies across the outline (the band's place, and its check)
        save_table(location.table(result.stack.pixel_size), out / f"{LOCATION_TABLE}.tsv")
        save_table(location.profiles(result.stack.pixel_size), out / f"{LOCATION_TABLE} - profiles.tsv")
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


def subtype_marker_names(settings: Settings, stack: EventStack) -> list[str]:
    """The names of the subtype markers as the tables write them ("<name>+"; settings resolved to numbers)."""
    return [stack.channel_names[stack.numbers.index(int(n))] for n in settings.channels.functional]


def start(settings: Settings) -> tuple[EventStack, Path]:
    """Load the input, make the output folder, start its log and save the settings."""
    stack = load(settings)
    moi_names = [stack.channel_names[stack.numbers.index(int(n))] for n in settings.channels.marker_of_interest]
    problems = selector_problems(settings.membrane, subtype_marker_names(settings, stack), moi_names)
    if problems:  # before anything is segmented
        raise ValueError("; ".join(problems))
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
                "doublets": len(result.doublets),
                "markers of interest": [
                    {"marker": m.name, "channel": m.number, "marker on": m.carried_by,
                     "medians (free membrane or whole cell)": m.medians,
                     "membrane pair": f"{' / '.join(m.cells)} with {' / '.join(m.partners)}",
                     "marker check": m.expression, "membranes": _membrane_info(m.membranes),
                     "file name suffix": m.suffix}
                    for m in result.markers],
                "cellpose": cells.attrs.get("cellpose", {}), "seconds": round(seconds, 1)}
    (result.out_dir / "run.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    log.info("Done in %.1f s: %d cells, %d single tumor + single immune cell events", run_info["seconds"],
             len(cells), len(result.doublets))


def _membrane_info(mem: MembraneResult | None) -> dict[str, Any] | None:
    return None if mem is None else {"cells": len(mem.table), "with_interface": int(mem.has_interface.sum()),
                                     "control_window_px": mem.window}


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
