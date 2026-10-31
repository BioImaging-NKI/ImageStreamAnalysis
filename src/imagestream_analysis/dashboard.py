"""The dashboard of one dataset: `dashboard.html` with the results embedded, and the event images as image sheets in
`dashboard_files/` (loaded with <img>, so the page also works opened from disk, file://).

Event images: per channel, 8-bit, display range from the channel's background (median) to its 99.8th percentile,
as sheets of SHEET_COLUMNS x SHEET_ROWS events (WebP). Cells: outlines (polygons) per cell; interfaces: outlines
per event. The page composes the channels itself (colour per channel, on/off, brightness).
"""

from __future__ import annotations

import base64
import json
import math
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import measure

from . import __version__
from .config import CellOfInterest, Settings

if TYPE_CHECKING:
    from .pipeline import Result

SHEET_COLUMNS, SHEET_ROWS = 10, 10
FILES = "dashboard_files"
ROLE_LUTS = {"brightfield": "gray", "immune": "azure", "tumor": "green", "marker_of_interest": "magenta"}
OTHER_LUTS = ["yellow", "red", "cyan", "orange"]


def _clean(v: Any) -> Any:
    """JSON-safe: NaN/inf -> None, numpy -> Python."""
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    if isinstance(v, (float, np.floating)):
        f = float(v)
        return f if math.isfinite(f) else None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def _column(values: Any, digits: int = 4) -> list[Any]:
    arr = np.asarray(values)
    if arr.dtype.kind in "fc":
        return [None if not np.isfinite(x) else round(float(x), digits) for x in arr]
    if arr.dtype.kind in "iub":
        return [int(x) for x in arr]
    return [str(x) for x in arr]


def float32_b64(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a, dtype="<f4").tobytes()).decode("ascii")


def _outline(mask: np.ndarray, y0: int = 0, x0: int = 0) -> list[int]:
    """The outer boundary of a mask (a crop at y0, x0 of the event) as a flat [x0, y0, x1, y1, ...] list in half
    pixels of the event (pixel edges are at integer coordinates in the page)."""
    contours = measure.find_contours(np.pad(mask, 1).astype(np.float32), 0.5)  # type: ignore[no-untyped-call]
    if not contours:
        return []
    c = max(contours, key=len)
    # (row, col) of the padded crop, at pixel centres -> page coordinates (x, y) with pixel edges at integers,
    # in half pixels (marching squares on a mask only gives half-pixel positions)
    xy = np.round(np.column_stack([c[:, 1] - 0.5 + x0, c[:, 0] - 0.5 + y0]) * 2).astype(np.int64)
    return [int(v) for v in _drop_collinear(xy).ravel()]


def _drop_collinear(xy: np.ndarray) -> np.ndarray:
    """A closed polygon without repeated points and without points on a straight line between their neighbours."""
    xy = xy[np.any(xy != np.roll(xy, 1, axis=0), axis=1)]
    if len(xy) < 4:
        return xy
    before, after = xy - np.roll(xy, 1, axis=0), np.roll(xy, -1, axis=0) - xy
    turn = before[:, 0] * after[:, 1] - before[:, 1] * after[:, 0] != 0
    out: np.ndarray = xy[turn]
    return out


def cell_outlines(labels: np.ndarray, frames: np.ndarray, cell_labels: np.ndarray) -> list[list[int]]:
    """The outline of every cell, traced on its bounding box only."""
    boxes: dict[int, list[tuple[slice, ...] | None]] = {}
    out: list[list[int]] = []
    for frame, label in zip(frames, cell_labels, strict=True):
        if frame not in boxes:
            boxes[frame] = ndi.find_objects(labels[frame])
        sl = boxes[frame][label - 1] if label - 1 < len(boxes[frame]) else None
        if sl is None:
            out.append([])
            continue
        out.append(_outline(labels[frame][sl] == label, sl[0].start, sl[1].start))
    return out


def _outlines_of_masks(mask: np.ndarray) -> list[list[int]]:
    """All boundaries of a mask (several pieces), as flat half-pixel lists."""
    lab = measure.label(mask)  # type: ignore[no-untyped-call]
    out = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is not None and (o := _outline(lab[sl] == i, sl[0].start, sl[1].start)):
            out.append(o)
    return out


def _display_range(channel: np.ndarray, background: float, brightfield: bool) -> tuple[float, float]:
    sample = channel[:: max(1, len(channel) // 300)]  # a few hundred events are enough
    if brightfield:  # brightfield: around its background
        lo, hi = np.percentile(sample, [0.2, 99.8])
        return float(lo), float(hi)
    hi = float(np.percentile(sample, 99.8))
    return background, max(hi, background + 1)


def write_sheets(result: Result, settings: Settings, folder: Path) -> list[dict[str, Any]]:
    """The image sheets of all loaded channels; returns their description for the page."""
    stack = result.stack
    roles: dict[int, str] = {}  # a channel with two roles shows the first (marker of interest before functional)
    for role, refs in settings.channels.roles().items():
        for ref in refs:
            roles.setdefault(stack.number(ref), role)
    numbers = stack.numbers
    luts_used: set[str] = set()
    t, _, h, w = stack.data.shape
    per_sheet = SHEET_COLUMNS * SHEET_ROWS
    channels = []
    jobs: list[tuple[np.ndarray, Path]] = []
    others = iter(OTHER_LUTS * 4)
    for i, (number, name) in enumerate(zip(numbers, stack.channel_names, strict=True)):
        role = roles.get(number, "")
        data = stack.data[:, i].astype(np.float32)
        lo, hi = _display_range(data, stack.medians[i] if stack.medians else 0.0, role == "brightfield")
        scaled = np.clip((data - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)
        for s in range(math.ceil(t / per_sheet)):
            part = scaled[s * per_sheet : (s + 1) * per_sheet]
            sheet = np.zeros((SHEET_ROWS * h, SHEET_COLUMNS * w), dtype=np.uint8)
            for k, im in enumerate(part):
                r, c = divmod(k, SHEET_COLUMNS)
                sheet[r * h : (r + 1) * h, c * w : (c + 1) * w] = im
            rows_used = math.ceil(len(part) / SHEET_COLUMNS)
            jobs.append((sheet[: rows_used * h], folder / f"c{number}_{s}.webp"))
        lut = ROLE_LUTS.get(role, "")
        if not lut or lut in luts_used:  # a second immune marker gets a colour of its own
            lut = next(others)
        luts_used.add(lut)
        channels.append({"number": number, "name": name, "role": role,
                         "lut": lut, "on": role in ROLE_LUTS,
                         "gain": 0.5 if role == "brightfield" else 1.0,  # brightfield under the fluorescence
                         "range": [round(lo, 2), round(hi, 2)]})
    with ThreadPoolExecutor(min(8, os.cpu_count() or 4)) as pool:  # Pillow encodes outside the GIL
        list(pool.map(lambda job: Image.fromarray(job[0]).save(job[1], quality=82, method=1), jobs))
    return channels


def build_data(result: Result, settings: Settings, channels: list[dict[str, Any]]) -> dict[str, Any]:
    stack = result.stack
    t, _, h, w = stack.data.shape
    cells = result.cells
    labels = result.labels
    ids = stack.object_ids
    of_interest = "Tumor" if settings.channels.cell_of_interest == CellOfInterest.TUMOR_CELL else "Immune"

    outlines = cell_outlines(labels, cells["frame"].to_numpy(), cells["label"].to_numpy())
    interface = result.masks.immune_interface | result.masks.tumor_interface
    interfaces = [_outlines_of_masks(interface[e]) if interface[e].any() else [] for e in range(t)]

    ev = result.events
    data: dict[str, Any] = {
        "version": __version__,
        "name": stack.name,
        "pixel_size": stack.pixel_size,
        "tile": [w, h],
        "n_events": t,
        "n_events_in_file": stack.n_events_in_file,
        "sheet": [SHEET_COLUMNS, SHEET_ROWS],
        "files": FILES,
        "channels": channels,
        "cell_of_interest": of_interest,
        "events": {
            "object": _column(ids) if ids is not None else None,
            "t_cells": _column(ev["Immune cells"]), "tumor_cells": _column(ev["Tumor cells"]),
            "t_interface": _column(ev["mean immune interface"]), "t_rest": _column(ev["mean immune NOT in interface"]),
            "tu_interface": _column(ev["mean tumor interface"]), "tu_rest": _column(ev["mean tumor NOT in interface"]),
            "t_ratio": _column(ev["Immune membrane ratio"]), "tu_ratio": _column(ev["Tumor membrane ratio"]),
            "interfaces": interfaces,
        },
        "cells": {
            "event": _column(cells["frame"]), "label": _column(cells["label"]),
            "type": _column(cells["Cell type"]), "area": _column(cells["Area"]),
            "t_int": _column(cells["Immune intensity"]), "tu_int": _column(cells["Tumor intensity"]),
            "ratio": _column(cells["Marker ratio"]), "outline": outlines,
        },
        "membranes": None,
    }
    mem = result.membranes
    if mem is not None:
        tab = mem.table
        index = {(int(f) - 1, int(lab)): i for i, (f, lab) in enumerate(zip(cells["frame"] + 1, cells["label"],
                                                                                strict=True))}
        data["membranes"] = {
            "cell": [index.get((int(f) - 1, int(lab)), -1) for f, lab in zip(tab["Frame"], tab["Label"],
                                                                               strict=True)],
            "sd": _column(tab["stdDevs above mean"]), "ratio": _column(tab["Ratio interface/non-interface"]),
            "control": _column(tab["control stdDevs above mean"]), "interface_length": _column(tab["interface length"]),
            "length": _column(mem.lengths), "width": int(mem.profiles.shape[1]),
            "profiles": float32_b64(mem.profiles),
            "average": _column(mem.average["profile"]),
            "hist_bin_start": _column(mem.histogram_interface["bin_start"]),
            "hist_interface": _column(mem.histogram_interface["counts"]),
            "hist_control": _column(mem.histogram_control["counts"]),
            "control_windows": len(mem.control_scores), "window": mem.window,
            "bin_width": settings.membrane.histogram_bin_width,
        }
    return data


def write_dashboard(result: Result, settings: Settings) -> Path:
    out = result.out_dir
    folder = out / FILES
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    channels = write_sheets(result, settings, folder)
    data = build_data(result, settings, channels)
    template = resources.files("imagestream_analysis").joinpath("templates", "dashboard.html").read_text("utf-8")
    payload = json.dumps(_clean(data), separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    path = out / "dashboard.html"
    path.write_text(template.replace("__ISA_VERSION__", __version__).replace('"__ISA_DATA__"', payload),
                    encoding="utf-8")
    return path
