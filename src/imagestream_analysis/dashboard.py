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
from .config import Settings
from .interfaces import moi_column

if TYPE_CHECKING:
    from .membrane import Location, Profiles
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
        if all(type(x) is int or type(x) is str for x in v):  # outlines, labels: nothing to clean (fast)
            return list(v)
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
    sample = channel[:: max(1, len(channel) // 300)].astype(np.float32)  # a few hundred events are enough
    if brightfield:  # brightfield: around its background
        lo, hi = np.percentile(sample, [0.2, 99.8])
        return float(lo), float(hi)
    hi = float(np.percentile(sample, 99.8))
    return background, max(hi, background + 1)


def _to_bytes(data: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """The display range lo..hi as 0..255. Integer images of up to 16 bits through a lookup table (2026-11-05: the
    same bytes as the float computation, without a float copy of the whole channel)."""
    if data.dtype.kind in "ui" and data.dtype.itemsize <= 2:
        offset = int(np.iinfo(data.dtype).min)
        values = np.arange(offset, offset + 2 ** (8 * data.dtype.itemsize), dtype=np.float32)
        lut = np.clip((values - np.float32(lo)) / np.float32(hi - lo) * np.float32(255.0), 0, 255).astype(np.uint8)
        out: np.ndarray = lut[data.astype(np.int64) - offset] if offset else lut[data]
        return out
    scaled: np.ndarray = np.clip((data.astype(np.float32) - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)
    return scaled


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
        raw = stack.data[:, i]
        lo, hi = _display_range(raw, stack.medians[i] if stack.medians else 0.0, role == "brightfield")
        scaled = _to_bytes(raw, lo, hi)
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

    outlines = cell_outlines(labels, cells["frame"].to_numpy(), cells["label"].to_numpy())
    interface = result.masks.immune_interface | result.masks.tumor_interface
    interfaces = [_outlines_of_masks(interface[e]) if interface[e].any() else [] for e in range(t)]

    ev = result.events
    row_of = {int(lab): r for r, lab in enumerate(cells["label"])}
    data: dict[str, Any] = {
        "version": __version__,
        "name": stack.name,
        "pixel_size": stack.pixel_size,
        "tile": [w, h],
        "n_events": t,
        "n_events_in_file": stack.n_events_in_file,
        "sources": stack.source_names or [],
        "sheet": [SHEET_COLUMNS, SHEET_ROWS],
        "files": FILES,
        "channels": channels,
        "has_moi": bool(result.markers),
        "artefact_min_positive": settings.classification.artefact_min_positive,
        "events": {
            "object": _column(ids) if ids is not None else None,
            "source": _column(stack.sources) if stack.sources is not None else None,  # the part of each event
            "t_cells": _column(ev["Immune cells"]), "tumor_cells": _column(ev["Tumor cells"]),
            "interface": _column(ev["Tumor-immune interface"].astype(int)),
            # the macro's per-event measure, for every marker of interest (interfaces.event_table)
            "moi": [{"name": mk.name, "carried_by": mk.carried_by or "Tumor",
                     **{key: _column(ev[moi_column(col, mk.name)])
                        for key, col in [("t_interface", "mean immune interface"),
                                         ("t_rest", "mean immune NOT in interface"),
                                         ("tu_interface", "mean tumor interface"),
                                         ("tu_rest", "mean tumor NOT in interface"),
                                         ("t_ratio", "Immune membrane ratio"), ("tu_ratio", "Tumor membrane ratio")]}}
                    for mk in result.markers],
            "interfaces": interfaces,
        },
        "cells": {
            "event": _column(cells["frame"]), "label": _column(cells["label"]),
            "type": _column(cells["Cell type"]), "subtype": _column(cells["Subtype"]), "area": _column(cells["Area"]),
            "t_int": _column(cells["Immune intensity"]), "tu_int": _column(cells["Tumor intensity"]),
            "ratio": _column(cells["Marker ratio"]), "outline": outlines,
        },
        "markers": markers_data(result, settings),
        "pairs": [row_of[int(v)] for v in result.pairs.ravel()],  # cells that touch, as rows of the cell table
        "membranes": None,
    }
    measured = [mk for mk in result.markers if mk.profiles is not None]
    if measured:  # every cell's membrane: the page computes the enrichment of any pair, for each marker of interest
        prof = measured[0].profiles
        assert prof is not None
        m = settings.membrane
        # the outlines and neighbour tracks are the same for every marker: once; a profile per marker
        track_rows = [np.array([row_of.get(int(v), -1) for v in t], dtype=np.int32) for t in prof.tracks]
        data["membranes"] = {
            "cell": [row_of[int(lab)] for lab in prof.labels],
            "length": [len(p) for p in prof.profiles],
            "track": base64.b64encode(np.concatenate(track_rows).astype("<i4").tobytes() if track_rows
                                      else b"").decode("ascii"),
            "line": membrane_lines(prof), "band_width": prof.width, "location": location_data(prof.location),
            "markers": [{"name": mk.name, "carried_by": mk.carried_by, "cells": mk.cells, "partners": mk.partners,
                         "expression": mk.expression,
                         "profile": float32_b64(np.concatenate(mk.profiles.profiles) if mk.profiles.profiles
                                                else np.zeros(0)),
                         # the sum across the band (the dashboard's switch; additive at an interface)
                         "profile_sum": float32_b64(np.concatenate(mk.profiles.sums) if mk.profiles.sums
                                                    else np.zeros(0))}
                        for mk in measured if mk.profiles is not None],
            "control_repeats": m.control_repeats, "max_histogram_value": m.max_histogram_value,
            "bin_width": m.histogram_bin_width, "inside_only": m.inside_only,
        }
    return data


LINE_SCALE = 8  # the straightening lines in 1/8 pixel


def location_data(loc: Location | None) -> dict[str, Any] | None:
    """Where the membrane lies across the outline (membrane.Location), for the membrane card: the average profile of
    each marker per cell type and place (pixels; + = outside the cell), and where the band was centred."""
    if loc is None:
        return None
    return {"distance": [round(float(d), 3) for d in loc.distance], "moved": loc.moved, "mode": loc.mode,
            "own": int(np.isfinite(loc.cell_free).sum()),
            "radial": [{"type": t, "at": w, "marker": m, "cells": loc.cells[(t, w, m)], "mean": _column(v, 3)}
                       for (t, w, m), v in sorted(loc.radial.items())],
            "offset": [{"type": t, "at": w, "offset": round(float(loc.offset.get((t, w), 0.0)), 3),
                        "source": loc.source[(t, w)]} for (t, w) in sorted(loc.source)]}


def membrane_lines(prof: Profiles) -> str:
    """The straightening line of every membrane (one point per profile sample, in the profiles' order) as base64 of
    little-endian int16 x, y pairs in 1/LINE_SCALE pixel of the event: the page draws the band that was sampled."""
    if not prof.lines:
        return ""
    xy = np.concatenate(prof.lines)[:, ::-1] * LINE_SCALE
    return base64.b64encode(np.round(xy).clip(-32768, 32767).astype("<i2").tobytes()).decode("ascii")


def markers_data(result: Result, settings: Settings) -> list[dict[str, Any]]:
    """The functional (subtype) markers: per cell the background-subtracted mean, and the thresholds of the run."""
    cells, th = result.cells, result.thresholds
    out = []
    for number in settings.channels.functional:
        prefix = f"Mean ch{result.stack.number(number)} ("
        column = next(c for c in cells.columns if c.startswith(prefix))
        name = column[len(prefix) : -1]
        rows = th[th["marker"] == name]
        thr = {str(c): float(v) for c, v in zip(rows["cells"], rows["threshold"], strict=True)}
        how = list(rows.loc[rows["cells"] == "immune cells", "how"]) if "how" in rows else []
        out.append({"name": name, "values": _column(cells[column], 3), "threshold": thr.get("immune cells"),
                    "threshold_all": thr.get("artefact rule"),
                    "method": str(rows["method"].iloc[0]) if len(rows) else "",
                    "how": str(how[0]) if how else ""})
    return out


def page_template(name: str) -> str:
    """A page of templates/ with the context help filled in: help.css, help.js and help.json (made from the manual by
    tools/manual/make_help.py), at the placeholders /*__ISA_HELP_CSS__*/, /*__ISA_HELP_JS__*/ and __ISA_HELP__."""
    def read(name: str) -> str:
        return resources.files("imagestream_analysis").joinpath("templates", name).read_text(encoding="utf-8")
    return (read(name).replace("/*__ISA_HELP_CSS__*/", read("help.css")).replace("/*__ISA_HELP_JS__*/", read("help.js"))
            .replace("__ISA_HELP__", read("help.json").replace("</", "<\\/")).replace("__ISA_LOGO__", logo_url()))


def logo_url() -> str:
    """The logo as a data URL (favicon and header of the pages; templates/logo.svg, made by tools/make_logo.py)."""
    svg = resources.files("imagestream_analysis").joinpath("templates", "logo.svg").read_bytes()
    return "data:image/svg+xml;base64," + base64.b64encode(svg).decode("ascii")


def write_dashboard(result: Result, settings: Settings) -> Path:
    out = result.out_dir
    folder = out / FILES
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    channels = write_sheets(result, settings, folder)
    data = build_data(result, settings, channels)
    template = page_template("dashboard.html")
    payload = json.dumps(_clean(data), separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    path = out / "dashboard.html"
    path.write_text(template.replace("__ISA_VERSION__", __version__).replace('"__ISA_DATA__"', payload),
                    encoding="utf-8")
    return path
