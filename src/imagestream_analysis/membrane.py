"""Membrane profiles and the enrichment at the interface (macro lines 700-996), between two cell types of choice.

Every cell (artefacts excepted) gets its membrane profile of each marker of interest and a neighbour track: for each
position along the outline, the label of the nearest other cell within the interface reach (0 = none). A pair of
cell types (`membrane.cells`, `membrane.partners`) then only chooses which positions are the interface: those that
touch a partner; the rest of the membrane is the part that touches no cell. (Decided 2026-10-31: with several cell
types the interfaces between all of them are kept by position, and any pair is a selection.) The outlines and tracks
are the same for every marker of interest; only the profiles differ.

Without a pair in the settings, the measured cells are the cell type (Tumor or Immune) that carries the marker, decided
per marker from the data (`carrier`; the user, 2026-11-03), against the other type.

Deviations from the macro, decided 2026-10-22 (SPEC I8-I11):
- one membrane per cell (the macro: one per blob of touching cells of the same class);
- the length of each profile is kept explicitly (the macro marks the end with NaN, which also hits true zeros);
- the control for cells without an interface: `control_repeats` random windows per cell (unseeded), as long as the
  median interface length (the macro: one window of 25 px);
- the average profile only uses cells with an interface.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field
from multiprocessing.shared_memory import SharedMemory
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure

from .config import MembraneSettings
from .interactions import type_label

log = logging.getLogger(__name__)
INTERFACE_OFFSET = 20  # the aligned interfaces are centred here (macro: interfaceOffset)
HISTOGRAM_MIN = -5.0  # macro: the histograms start at -5 standard deviations


@dataclass
class MembraneResult:
    table: pd.DataFrame  # one row per cell of interest
    profiles: np.ndarray  # (cells, max length) normalized MOI profiles, aligned on the interface; NaN = beyond the end
    lengths: np.ndarray  # length of each profile
    has_interface: np.ndarray  # bool per cell
    average: pd.DataFrame  # average aligned profile of the cells with an interface
    histogram_interface: pd.DataFrame
    histogram_control: pd.DataFrame
    control_scores: np.ndarray  # all control windows of all cells without interface (the bootstrap distribution)
    window: int  # control window length (pixels along the membrane)
    # the same profiles minus the partner's free membrane at the interface positions (excess over both membranes)
    excess_profiles: np.ndarray | None = None
    control_log2: np.ndarray | None = None  # log2 ratio of the same control windows (LOG2)


EXCESS = "excess stdDevs above both membranes"  # the column of the excess (partner_excess) in the membrane table
# the score minus the median score of the measured cells at contacts with cells of their own selection (2026-11-05):
# an empirical reference for two stained membranes side by side, where the excess's assumption failed
SAME = "stdDevs above same-type contacts"
MIN_SAME = 10  # fewer measured cells with a same-type contact: no reference (NaN)
# the main readout (2026-11-05): log2 of interface mean / mean of the own membrane that touches no cell (both
# background-subtracted). A fold change: independent of how patchy or noisy the free membrane is, and symmetric.
# NaN when either mean is not above 0 (a cell that hardly carries the marker).
LOG2 = "log2 ratio interface/non-interface"
LOG2_SAME = "log2 ratio above same-type contacts"  # LOG2 minus the median LOG2 at same-type contacts (as SAME)


def free_means(membranes: Profiles) -> dict[int, float]:
    """The mean of every cell's membrane that touches no cell (background-subtracted), by label; NaN without one.
    The same per-cell value as the marker check (per_cell_intensity)."""
    return {int(lab): float(p[t == 0].mean()) if (t == 0).any() else float("nan")
            for lab, p, t in zip(membranes.labels, membranes.profiles, membranes.tracks, strict=True)}


def partner_excess(p: np.ndarray, f: np.ndarray, track: np.ndarray, free_of: dict[int, float]
                   ) -> tuple[float, np.ndarray]:
    """The excess over both membranes of one cell (decided 2026-11-03): at an interface the band holds both
    membranes, closer than the resolution, and background-subtracted intensities add. So

        excess = (mean(p[interface]) - own free mean - mean over the interface positions of the partner's free mean)
                 / own free sd

    = the score (stdDevs above mean) minus the mean partner free intensity / own free sd. Returns that subtracted
    mean (in intensity, NaN without an interface or without any partner with a free membrane) and, per position,
    the partner's free mean to subtract (0 outside the interface; a partner without a free membrane: the mean of the
    other positions)."""
    sub = np.zeros(len(p))
    if not f.any():
        return float("nan"), sub
    vals = np.array([free_of.get(int(lab), np.nan) for lab in track[f]], dtype=float)
    if not np.isfinite(vals).any():
        return float("nan"), sub
    m = float(np.nanmean(vals))
    sub[f] = np.where(np.isfinite(vals), vals, m)
    return m, sub


def outline(mask: np.ndarray, smoothing: float = 2.0) -> np.ndarray | None:
    """The outer boundary of a cell as a closed line, resampled at 1-pixel steps: (n, 2) array of (y, x).

    Marching squares at level 0.5 (the edge between the cell's pixels and the background), smoothed along the line
    (Gaussian, sigma in samples) as Fiji smooths the traced outline with a spline before straightening."""
    contours = measure.find_contours(np.pad(mask, 1).astype(np.float32), 0.5)  # type: ignore[no-untyped-call]
    if not contours:
        return None
    c = max(contours, key=len)[:-1] - 1.0  # closed: last point = first; undo the padding
    if len(c) < 4:
        return None
    c = ndi.gaussian_filter1d(c, smoothing, axis=0, mode="wrap")
    seg = np.linalg.norm(np.diff(np.vstack([c, c[:1]]), axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    perimeter = s[-1]
    n = max(int(round(perimeter)), 4)
    t = np.arange(n) * perimeter / n
    closed = np.vstack([c, c[:1]])
    return np.column_stack([np.interp(t, s, closed[:, 0]), np.interp(t, s, closed[:, 1])])


def straighten(image: np.ndarray, line: np.ndarray, width: int, order: Literal[0, 1] = 1) -> np.ndarray:
    """Fiji Straighten along a closed line: `width` samples across the line (bilinear; order 0: nearest), (width,
    n)."""
    tangent = (np.roll(line, -1, axis=0) - np.roll(line, 1, axis=0)) / 2
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-9)
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    offsets = np.arange(width) - (width - 1) / 2
    pts = line[None] + offsets[:, None, None] * normal[None]
    return ndi.map_coordinates(image, [pts[..., 0], pts[..., 1]], order=order, mode="nearest")


def _runs(flags: np.ndarray) -> list[tuple[int, int]]:
    """Runs of True in a circular array: (start, length), the run over the end counted once."""
    n = len(flags)
    if flags.all():
        return [(0, n)]
    if not flags.any():
        return []
    start = int(np.argmin(flags))  # a False position: no run wraps around it
    rolled = np.roll(flags, -start)
    edges = np.diff(np.concatenate([[0], rolled.astype(np.int8), [0]]))
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    return [((int(a) + start) % n, int(b - a)) for a, b in zip(starts, ends, strict=True)]


def _log2(window: np.ndarray, rest: np.ndarray) -> float:
    """log2(window mean / rest mean); NaN when either is empty or not above 0 (LOG2)."""
    if len(window) == 0 or len(rest) == 0:
        return float("nan")
    a, b = float(window.mean()), float(rest.mean())
    return float(np.log2(a / b)) if a > 0 and b > 0 else float("nan")


def _score(window: np.ndarray, rest: np.ndarray) -> float:
    """(window mean - rest mean) / rest sd (sample sd, as ImageJ)."""
    if len(window) == 0 or len(rest) < 2:
        return float("nan")
    sd = rest.std(ddof=1)
    return float((window.mean() - rest.mean()) / sd) if sd > 0 else float("nan")


HISTOGRAM_PERCENTILE = 99.0  # automatic maximum: this percentile of the interface or the control values
HISTOGRAM_MIN_TOP = 5.0  # automatic maximum: at least this
HISTOGRAM_MAX_BINS = 60  # automatic bin width: the smallest nice step that gives at most this many bins
HISTOGRAM_FALLBACK_TOP = 20.0  # automatic maximum when there are no values (the macro's default)


def _nice_step(raw: float) -> float:
    """The smallest of 0.1, 0.2, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, ... that is at least `raw` (each divides 5)."""
    if raw <= 0.1:
        return 0.1
    for step in (0.2, 0.25, 0.5):
        if raw <= step + 1e-12:
            return step
    unit = 1.0
    while True:
        for m in (1.0, 2.5, 5.0):
            if raw <= m * unit + 1e-12:
                return m * unit
        unit *= 10


def histogram_scale(interface: np.ndarray, control: np.ndarray, max_value: float | None = None,
                    bin_width: float | None = None) -> tuple[float, float, int]:
    """The histogram range [-5, top], its bin width and number of bins. An empty (None or 0) setting is automatic:

    - top: the 99th percentile (numpy's linear interpolation) of the finite interface values or of the control
      windows, whichever is higher, at least 5 (20 without values), rounded up to a whole number of bins;
    - bin width: the smallest nice step (_nice_step) that gives at most 60 bins over [-5, top].

    Explicit values are used as given, as in the macro: round((top + 5) / bin width) bins. The dashboard repeats this
    rule in JS (histScale in templates/dashboard.html)."""
    if max_value:
        top = float(max_value)
    else:
        finite = [v[np.isfinite(v)] for v in (np.asarray(interface, float), np.asarray(control, float))]
        pcts = [float(np.percentile(v, HISTOGRAM_PERCENTILE)) for v in finite if len(v)]
        top = max(HISTOGRAM_MIN_TOP, *pcts) if pcts else HISTOGRAM_FALLBACK_TOP
    width = float(bin_width) if bin_width else _nice_step((top - HISTOGRAM_MIN) / HISTOGRAM_MAX_BINS)
    if not max_value:
        top = round(HISTOGRAM_MIN + float(np.ceil((top - HISTOGRAM_MIN) / width - 1e-9)) * width, 10)
    n_bins = max(1, int(round((top - HISTOGRAM_MIN) / width)))
    return top, width, n_bins


def _histogram(values: np.ndarray, top: float, n_bins: int) -> pd.DataFrame:
    """As the macro's Distribution: n_bins bins over [-5, top]; values outside are not counted."""
    counts, edges = np.histogram(values[np.isfinite(values)], bins=n_bins, range=(HISTOGRAM_MIN, top))
    return pd.DataFrame({"bin_start": edges[:-1], "counts": counts})


@dataclass
class Location:
    """Where the membrane lies relative to the outline (2026-11-04): per cell type of the measured cell and where
    (free membrane: "", or a contact with a cell type), the marker of interest averaged across the outline over the
    cells (`radial`, at the distances `distance` in pixels, + = outside the cell), its peak (`offset`, pixels), the
    marker it was taken from, and the spread of the peaks of single cells (`cell_offsets`)."""

    distance: np.ndarray
    radial: dict[tuple[str, str, str], np.ndarray] = field(default_factory=dict)  # (type, where, marker) -> mean
    cells: dict[tuple[str, str, str], int] = field(default_factory=dict)
    offset: dict[tuple[str, str], float] = field(default_factory=dict)  # (type, where) -> pixels, as used
    source: dict[tuple[str, str], str] = field(default_factory=dict)  # the marker it comes from, or why none
    cell_offsets: dict[tuple[str, str], np.ndarray] = field(default_factory=dict)  # per cell, from that marker
    # per cell (index in the measured cells), the peak of its own free membrane (pixels) from the marker that placed
    # its type's free band; NaN: too little free membrane, or no peak of its own (centre "cell" then uses the type's)
    cell_free: np.ndarray = field(default_factory=lambda: np.zeros(0))
    cell_types: list[str] = field(default_factory=list)  # the cell type of each measured cell (index as cell_free)
    mode: str = "outline"  # where the band was centred: "outline", "membrane" (per cell type) or "cell"

    @property
    def moved(self) -> bool:
        return self.mode != "outline"

    def shift(self, index: int, typ: str, where: np.ndarray) -> np.ndarray:
        """The distance (pixels, + = outward) from the outline to the band's centre at every position of one cell.
        "membrane": the offset of its type at each place (free, or a contact with a cell type; 0 where none was
        found). "cell": at the free membrane the cell's own free peak; at a contact that peak plus its type's
        difference between that contact and the free membrane (contacts are too short, and too close to what is
        tested, to place a band per cell); a cell without a peak of its own as "membrane"."""
        if not self.moved:
            return np.zeros(len(where))
        base = [self.offset.get((typ, str(w)), 0.0) for w in where]
        own = self.cell_free[index] if self.mode == "cell" and index < len(self.cell_free) else np.nan
        free = self.offset.get((typ, ""))
        if not np.isfinite(own) or free is None:
            return np.array(base)
        return np.array([own + (b - free if (typ, str(w)) in self.offset else 0.0) if w else own
                         for w, b in zip(where, base, strict=True)])

    def table(self, pixel_size: float) -> pd.DataFrame:
        """One row per cell type, place and marker: the peak of the average profile across the outline (µm, + =
        outside), the cells, whether it placed the band, and the spread of the single cells' peaks."""
        rows = []
        for (typ, where, marker), prof in sorted(self.radial.items()):
            used = self.source.get((typ, where)) == marker
            spread = self.cell_offsets.get((typ, where)) if used else None
            fin = spread[np.isfinite(spread)] if spread is not None else np.zeros(0)
            q = np.percentile(fin, [25, 50, 75]) * pixel_size if len(fin) else [np.nan] * 3
            rows.append({"cell type": typ, "at": f"contact with {where}" if where else "free membrane",
                         "marker": marker, "cells": self.cells[(typ, where, marker)],
                         "peak (µm from the outline)": _peak(self.distance, prof) * pixel_size,
                         "band placed by it": used and self.moved,
                         "band at (µm)": (self.offset.get((typ, where), 0.0) if self.moved else 0.0) * pixel_size,
                         "band per cell": self.mode == "cell" and used,
                         "single cells p25 (µm)": q[0], "single cells median (µm)": q[1],
                         "single cells p75 (µm)": q[2]})
        return pd.DataFrame(rows)

    def height(self, typ: str, marker: str) -> float:
        """How bright the free membrane of a cell type is in a marker: the peak of its average profile across the
        outline above the ends of that profile (0 without a peak; NaN without cells). It does not depend on the
        band (its width or place): the carrier of a marker is decided with it (2026-11-05)."""
        prof = self.radial.get((typ, "", marker))
        if prof is None:
            return float("nan")
        if not np.isfinite(_peak(self.distance, prof)):
            return 0.0
        k = max(1, len(prof) // 8)
        return max(0.0, float(prof.max() - (prof[:k].mean() + prof[-k:].mean()) / 2))

    def carrier(self, marker: str) -> tuple[str, dict[str, float]]:
        """The cell type (Tumor or Immune) whose free membrane is the brightest in this marker, and the heights."""
        heights = {t: self.height(t, marker) for t in ("Tumor", "Immune")}
        finite = {t: h for t, h in heights.items() if np.isfinite(h)}
        return (max(finite, key=lambda t: finite[t]) if finite and max(finite.values()) > 0 else ""), heights

    def profiles(self, pixel_size: float) -> pd.DataFrame:
        """The average profiles across the outline, long format (for plots)."""
        return pd.DataFrame([{"cell type": t, "at": f"contact with {w}" if w else "free membrane", "marker": m,
                              "distance (µm)": round(float(d) * pixel_size, 3), "mean": float(v)}
                             for (t, w, m), prof in sorted(self.radial.items()) for d, v in zip(self.distance, prof,
                                                                                                   strict=True)])


LOCATE_UM = 3.0  # how far on each side of the outline the membrane is looked for
LOCATE_STEP = 0.25  # pixels between the samples across the outline when locating
MIN_LOCATE_CELLS = 20  # fewer cells (with at least MIN_SAMPLES positions there): the band is not moved there
MIN_SAMPLES = 5


def _peak(distance: np.ndarray, profile: np.ndarray) -> float:
    """The position of the maximum, refined with a parabola; NaN at the end of the window (no peak inside it)."""
    if not len(profile) or not np.isfinite(profile).all():
        return float("nan")
    k = int(np.argmax(profile))
    if k == 0 or k == len(profile) - 1:
        return float("nan")
    a, b, c = profile[k - 1], profile[k], profile[k + 1]
    curve = a - 2 * b + c
    step = float(distance[1] - distance[0])
    return float(distance[k] + (0.5 * (a - c) / curve * step if curve < 0 else 0.0))


def _prominence(profile: np.ndarray) -> float:
    """How clearly the average profile has a membrane: its peak above its ends, relative to its range (0-1)."""
    ends = (profile[: max(1, len(profile) // 8)].mean() + profile[-max(1, len(profile) // 8) :].mean()) / 2
    span = float(profile.max() - profile.min())
    return float((profile.max() - ends) / span) if span > 0 else 0.0


@dataclass
class _Cell:
    t: int
    label: int
    corner: tuple[int, int]
    box: tuple[slice, slice]
    own: np.ndarray
    line: np.ndarray  # the outline, (n, 2) y, x in the box
    outward: np.ndarray  # unit normals pointing out of the cell, (n, 2)
    track: np.ndarray  # nearest other cell within the reach (label, 0 = none)
    where: np.ndarray  # per position: "" (free) or the cell type of that neighbour
    typ: str  # the cell type of this cell


def _outward(own: np.ndarray, line: np.ndarray) -> np.ndarray:
    tangent = (np.roll(line, -1, axis=0) - np.roll(line, 1, axis=0)) / 2
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-9)
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    inside = ndi.map_coordinates(own.astype(np.float64), (line + 1.5 * normal).T, order=0, mode="constant")
    out: np.ndarray = -normal if inside.mean() > 0.5 else normal
    return out


def location_notes(loc: Location, width: int, pixel_size: float) -> list[tuple[int, str]]:
    """Log lines (level, text) on where the membrane lies: per cell type and place, the peak relative to the
    outline, the spread of single cells, and whether the band reaches it. A band on the outline (centre "outline")
    whose edge is less than a pixel beyond a membrane is a warning: at a free membrane it then reaches less of the
    membrane than at an interface, which fakes enrichment."""
    reach = (width - 1) / 2  # the outermost sample of the band, from its centre
    out: list[tuple[int, str]] = []
    for (typ, where), offset in sorted(loc.offset.items()):
        place = f"contacts with {where.lower()} cells" if where else "free membrane"
        cells = loc.cell_offsets.get((typ, where), np.zeros(0))
        cells = cells[np.isfinite(cells)]
        spread = (f"; single cells {np.percentile(cells, 25) * pixel_size:+.2f} to "
                  f"{np.percentile(cells, 75) * pixel_size:+.2f} µm (middle half)") if len(cells) else ""
        text = (f"Membrane of {typ.lower()} cells, {place}: {offset * pixel_size:+.2f} µm from the outline "
                f"(+ = outside; {loc.source[(typ, where)]}){spread}")
        if loc.mode == "cell":
            out.append((logging.INFO, text + ("; the band is centred on each cell's own free membrane ("
                                              f"{_own_count(loc, typ)})" if not where else
                                              "; at these contacts each cell's free place plus the type's difference")))
        elif loc.moved:
            out.append((logging.INFO, text + ": the band is centred there"))
        elif abs(offset) > reach - 1:
            out.append((logging.WARNING, text + f"; the band on the outline reaches {reach * pixel_size:.2f} µm on "
                        "each side and misses part of it: an interface may look enriched without cause "
                        "(membrane.band_centre = 'membrane', or a wider membrane.straightening_width)"))
        else:
            out.append((logging.INFO, text))
    said: set[str] = set()
    for (typ, where), why in sorted(loc.source.items()):
        if (typ, where) in loc.offset or where == "Artefact" or typ in said:
            continue
        if why.startswith("no marker of interest is carried"):  # one line for the whole cell type
            said.add(typ)
            out.append((logging.INFO, f"Membrane of {typ.lower()} cells: {why}"))
        else:
            out.append((logging.INFO, f"Membrane of {typ.lower()} cells, "
                        f"{f'contacts with {where.lower()} cells' if where else 'free membrane'}: {why}"))
    return out


CARRY_FRACTION = 0.5  # a cell type carries a marker at this fraction of the brightest type's median, or more
MIN_FREE = 20  # positions of free membrane a cell needs for a peak of its own (centre "cell")


def _own_count(loc: Location, typ: str) -> str:
    own = loc.cell_free[[i for i, t in enumerate(loc.cell_types) if t == typ]] if len(loc.cell_free) else np.zeros(0)
    return f"{int(np.isfinite(own).sum())} of {len(own)} cells by their own peak, the others by the type's"


def _relative(level: dict[str, float] | None, typ: str) -> float:
    """The marker's level in this cell type relative to the brightest type (1 without levels; 0 without a value)."""
    if level is None:
        return 1.0
    values = [v for v in level.values() if np.isfinite(v)]
    top = max(values) if values else float("nan")
    own = level.get(typ, float("nan"))
    return float(own / top) if np.isfinite(own) and np.isfinite(top) and top > 0 else 0.0


def _cells_of(labels: np.ndarray, t0: int, rows: Sequence[tuple[int, int]], type_of: dict[tuple[int, int], str],
              width: int, reach: int, half: float) -> list[_Cell]:
    """The geometry of these cells (rows: event, label; events t0.. of `labels`): outline, outward normals, the
    neighbour track (on the outline) and what each position touches."""
    pad = width + reach + int(np.ceil(half)) + 1  # room for the band, the neighbours and the search
    objects: dict[int, list[Any]] = {}
    out = []
    for t, label in rows:
        frame = labels[t - t0]
        if t not in objects:
            objects[t] = ndi.find_objects(frame)
        sl = objects[t][label - 1] if label - 1 < len(objects[t]) else None
        if sl is None:
            continue
        y0, x0 = max(sl[0].start - pad, 0), max(sl[1].start - pad, 0)
        box = (slice(y0, sl[0].stop + pad), slice(x0, sl[1].stop + pad))
        sub = frame[box]
        own = sub == label
        line = outline(own)
        if line is None:
            continue
        # every own pixel within `reach` of another cell: the label of the nearest one (as imagej.dilate: taxicab)
        others = np.where(own, 0, sub)
        dist, (iy, ix) = ndi.distance_transform_cdt(others == 0, metric="taxicab", return_indices=True)
        nearest = np.where(own & (dist <= reach), others[iy, ix], 0)
        track = straighten(nearest.astype(np.float64), line, width, order=0).max(axis=0).astype(np.int64)
        where = np.array([type_of.get((t, int(v)), "") if v else "" for v in track], dtype=object)
        out.append(_Cell(t, label, (y0, x0), box, own, line, _outward(own, line), track, where, type_of[(t, label)]))
    return out


Across = dict[tuple[str, str], np.ndarray]  # (place, marker) -> the marker averaged across the outline there


def _across(cells: list[_Cell], mois: Sequence[np.ndarray], t0: int, backgrounds: Sequence[float],
            names: Sequence[str], distance: np.ndarray) -> tuple[list[Across], list[dict[str, np.ndarray | None]]]:
    """Per cell, every marker sampled across the outline (`distance`, + outward) and averaged along each place with
    at least MIN_SAMPLES positions; and along its free membrane if it has MIN_FREE positions (else None)."""
    places: list[Across] = []
    frees: list[dict[str, np.ndarray | None]] = []
    for c in cells:
        pts = c.line[None] + distance[:, None, None] * c.outward[None]
        mine: Across = {}
        free_of: dict[str, np.ndarray | None] = {}
        for moi, background, name in zip(mois, backgrounds, names, strict=True):
            across = ndi.map_coordinates(moi[c.t - t0][c.box].astype(np.float64) - background,
                                         [pts[..., 0], pts[..., 1]], order=1, mode="nearest")
            for where in np.unique(c.where):
                at = c.where == where
                if at.sum() >= MIN_SAMPLES:
                    mine[(str(where), name)] = across[:, at].mean(axis=1)
            free = c.where == ""
            free_of[name] = across[:, free].mean(axis=1) if free.sum() >= MIN_FREE else None
        places.append(mine)
        frees.append(free_of)
    return places, frees


def locate(cells: list[_Cell], places: list[Across], frees: list[dict[str, np.ndarray | None]],
           names: Sequence[str], distance: np.ndarray, levels: Sequence[dict[str, float]] | None = None) -> Location:
    """From every cell's marker averaged across the outline (`_across`): the average per cell type and place (free,
    or a contact with a cell type), and where the band goes: for each type and place, the peak of the marker with
    the clearest membrane there (prominence), among the markers that this cell type carries: `levels` per marker,
    the median intensity per cell type; a type carries a marker at CARRY_FRACTION of the brightest type or more
    (None: all), and the marker most specific to the type (its level relative to the brightest) comes first, then
    the prominence. A marker that a cell type does not carry shows no membrane of it, or signal elsewhere in the
    cell (seen on real data, 2026-11-04). None clear, or too few cells: on the outline (0)."""
    loc = Location(distance)
    sums: dict[tuple[str, str, str], list[np.ndarray]] = {}
    for c, mine in zip(cells, places, strict=True):
        for (where, name), prof in mine.items():
            sums.setdefault((c.typ, where, name), []).append(prof)
    loc.cell_types = [c.typ for c in cells]
    loc.cell_free = np.full(len(cells), np.nan)
    for key, profs in sums.items():
        loc.radial[key] = np.mean(profs, axis=0)
        loc.cells[key] = len(profs)
    if levels is None:  # how bright each type's free membrane is in each marker (Location.height)
        types = sorted({t for t, _, _ in sums})
        levels = [{t: loc.height(t, m) for t in types} for m in names]
    for typ, where in sorted({(t, w) for t, w, _ in sums}):
        rel = {m: _relative(levels[k], typ) for k, m in enumerate(names)}
        own = [m for m in names if rel[m] >= CARRY_FRACTION]
        if not own:
            loc.source[(typ, where)] = "no marker of interest is carried by these cells: on the outline"
            continue
        candidates = [(rel[m], _prominence(loc.radial[(typ, where, m)]), m) for m in own
                      if (typ, where, m) in loc.radial and loc.cells[(typ, where, m)] >= MIN_LOCATE_CELLS
                      and np.isfinite(_peak(distance, loc.radial[(typ, where, m)]))]
        if not candidates:
            loc.source[(typ, where)] = "too few cells, or no peak: on the outline"
            continue
        candidates = [c for c in candidates if c[1] >= 0.3]
        if not candidates:
            loc.source[(typ, where)] = "no clear membrane: on the outline"
            continue
        marker = max(candidates)[2]  # the most specific marker of this type, then the clearest peak
        loc.offset[(typ, where)] = _peak(distance, loc.radial[(typ, where, marker)])
        loc.source[(typ, where)] = marker
        loc.cell_offsets[(typ, where)] = np.array([_peak(distance, p) for p in sums[(typ, where, marker)]])
        if not where:  # each cell's own free peak, from the same marker (centre "cell")
            for i, free_of in enumerate(frees):
                own_free = free_of.get(marker)
                if loc.cell_types[i] == typ and own_free is not None and _prominence(own_free) >= 0.3:
                    loc.cell_free[i] = _peak(distance, own_free)
    return loc


@dataclass
class Profiles:
    """The membrane of every cell: its profile (marker of interest minus background, maximum across the membrane)
    and its neighbour track (label of the nearest other cell within the reach, 0 = none), same length.

    sums: the same profiles as the sum across the band instead of the maximum (2026-11-03). At an interface both
    membranes lie in the band, and a sum adds them where a maximum does not: the projection for the excess over both
    membranes (a switch in the dashboard; the tables keep the maximum)."""

    labels: np.ndarray  # the cell (label) of each profile
    profiles: list[np.ndarray]
    tracks: list[np.ndarray]
    sums: list[np.ndarray] = field(default_factory=list)
    # the line along which each membrane was straightened, in pixels of its event ((n, 2) of y, x; one point per
    # profile sample), and the width of the band across it: the dashboard draws the band that was sampled
    lines: list[np.ndarray] = field(default_factory=list)
    width: int = 0
    location: Location | None = None  # where the membrane was found (centre "membrane"), shared by the markers
    _free: pd.Series | None = field(default=None, init=False, repr=False, compare=False)

    def free(self) -> pd.Series:
        """Every cell's mean of the membrane that touches no cell (NaN without one), by label; computed once (the
        marker check asks for it for every pair)."""
        if self._free is None:
            self._free = pd.Series([float(p[t == 0].mean()) if (t == 0).any() else np.nan
                                    for p, t in zip(self.profiles, self.tracks, strict=True)], index=self.labels)
        return self._free

    def summed(self) -> Profiles:
        """The same membranes with the sum across the band as profile."""
        return Profiles(self.labels, self.sums, self.tracks, self.sums, self.lines, self.width, self.location)


def membrane_profiles(labels: np.ndarray, cells: pd.DataFrame, moi: np.ndarray, background: float, width: int,
                      reach: int, inside_only: bool = False) -> Profiles:
    """Profiles and neighbour tracks of all cells but the artefacts, for one marker (membrane_profiles_of)."""
    return membrane_profiles_of(labels, cells, [moi], [background], width, reach, inside_only)[0]  # on the outline


def _bands(cells: list[_Cell], lines: list[np.ndarray], mois: Sequence[np.ndarray], t0: int,
           backgrounds: Sequence[float], width: int, inside_only: bool
           ) -> tuple[list[list[np.ndarray]], list[list[np.ndarray]]]:
    """Per marker, every cell's band straightened along its line: the maximum and the sum across it."""
    maxima: list[list[np.ndarray]] = [[] for _ in mois]
    sums: list[list[np.ndarray]] = [[] for _ in mois]
    for c, line in zip(cells, lines, strict=True):
        if inside_only:
            inside = straighten(c.own.astype(np.float64), line, width, order=0) > 0
            inside[(width - 1) // 2] |= ~inside.any(axis=0)  # a sliver of a cell: at least the middle sample
        for k, (moi, background) in enumerate(zip(mois, backgrounds, strict=True)):
            band = straighten(moi[c.t - t0][c.box].astype(np.float64) - background, line, width)
            if inside_only:
                band = np.where(inside, band, -np.inf)
            maxima[k].append(band.max(axis=0))  # macro: maximum projection across the membrane
            sums[k].append(np.where(np.isfinite(band), band, 0.0).sum(axis=0))  # additive: both membranes count
    return maxima, sums


PARALLEL_MIN_CELLS = 2000  # fewer cells: one process (starting workers costs more than it saves)


def _chunks(frames: np.ndarray, n: int) -> list[tuple[int, int]]:
    """Contiguous ranges of events [t0, t1) with about equal numbers of cells (frames: the event of every cell)."""
    if not len(frames):
        return []
    cuts = np.unique(np.quantile(frames, np.linspace(0, 1, n + 1)[1:-1]).astype(int))
    edges = [int(frames.min()), *[int(c) for c in cuts if frames.min() < c <= frames.max()], int(frames.max()) + 1]
    return [(a, b) for a, b in zip(edges[:-1], edges[1:], strict=True) if b > a]


@dataclass(frozen=True)
class _Shared:
    """An array in shared memory, by name: the worker processes read the stack there instead of receiving a pickled
    copy of their range of events."""

    name: str
    shape: tuple[int, ...]
    dtype: str


_ATTACHED: dict[str, tuple[SharedMemory, np.ndarray]] = {}  # in a worker: the arrays it attached (until it exits)


def _array(a: np.ndarray | _Shared) -> np.ndarray:
    if isinstance(a, np.ndarray):
        return a
    if a.name not in _ATTACHED:
        shm = SharedMemory(name=a.name)
        _ATTACHED[a.name] = (shm, np.ndarray(a.shape, dtype=a.dtype, buffer=shm.buf))
    return _ATTACHED[a.name][1]


def _share(a: np.ndarray) -> tuple[SharedMemory, _Shared]:
    shm = SharedMemory(create=True, size=max(a.nbytes, 1))
    np.ndarray(a.shape, dtype=a.dtype, buffer=shm.buf)[...] = a
    return shm, _Shared(shm.name, a.shape, a.dtype.str)


def _phase1(labels: np.ndarray | _Shared, mois: list[np.ndarray | _Shared], t0: int, t1: int,
            rows: list[tuple[int, int]], type_of: dict[tuple[int, int], str], width: int, reach: int, half: float,
            backgrounds: list[float], names: list[str], distance: np.ndarray
            ) -> tuple[list[_Cell], list[Across], list[dict[str, Any]]]:
    images = [_array(m)[t0:t1] for m in mois]
    cells = _cells_of(_array(labels)[t0:t1], t0, rows, type_of, width, reach, half)
    places, frees = _across(cells, images, t0, backgrounds, names, distance)
    return cells, places, frees


def _phase2(cells: list[_Cell], lines: list[np.ndarray], mois: list[np.ndarray | _Shared], t0: int, t1: int,
            backgrounds: list[float], width: int, inside_only: bool
            ) -> tuple[list[list[np.ndarray]], list[list[np.ndarray]]]:
    return _bands(cells, lines, [_array(m)[t0:t1] for m in mois], t0, backgrounds, width, inside_only)


def membrane_profiles_of(labels: np.ndarray, cells: pd.DataFrame, mois: Sequence[np.ndarray],
                         backgrounds: Sequence[float], width: int, reach: int,
                         inside_only: bool = False, centre: str = "outline", names: Sequence[str] | None = None,
                         pixel_size: float = 1.0, levels: Sequence[dict[str, float]] | None = None,
                         workers: int | None = None) -> list[Profiles]:
    """Profiles and neighbour tracks of all cells but the artefacts, one Profiles per marker (they share the labels
    and the tracks: the outlines are traced once). labels (T, Y, X); mois: (T, Y, X) raw markers of interest, their
    `backgrounds` are subtracted; reach: the interface reach (membrane.thickness), in 4-connected steps.

    The straightened band is `width` samples across a line along the cell. centre "outline" (macro 1.8): the line is
    the outline, half of the band outside the cell. centre "membrane" (2026-11-04): the line is moved to where the
    membrane lies (`locate`, per cell type, at the free membrane and at contacts with each cell type), because the
    outline lies outside a free membrane but on the contact line between touching cells; a band on the outline then
    reaches less of the free membrane than of an interface. centre "cell": each cell's free membrane on its own peak
    (`Location.shift`). The neighbour track is always taken on the outline.
    inside_only: the maximum is taken over the samples inside the cell only, so that the membrane of a partner that
    also carries the marker is not measured as this cell's.

    The cells are measured in worker processes, a range of events each (2026-11-04: this loop was the slowest part
    of a run), from PARALLEL_MIN_CELLS cells on; `workers` (default: the cores, at most 8). The results are the same
    as in one process, in the same order."""
    names = list(names) if names is not None else [f"marker {k + 1}" for k in range(len(mois))]
    roi = cells[cells["Cell type"] != "Artefact"]
    rows = list(zip(roi["frame"].astype(int), roi["label"].astype(int), strict=True))
    type_of = {(int(t), int(lab)): str(c) for t, lab, c in zip(cells["frame"], cells["label"], cells["Cell type"],
                                                               strict=True)}
    half = LOCATE_UM / pixel_size  # the search for the membrane, in pixels on each side of the outline
    distance = np.arange(-half, half + LOCATE_STEP / 2, LOCATE_STEP)
    n_workers = workers if workers is not None else min(os.cpu_count() or 1, 8)
    if len(rows) < PARALLEL_MIN_CELLS:
        n_workers = 1
    frames = np.array([t for t, _ in rows], dtype=np.int64)
    ranges = _chunks(frames, n_workers * 3) if n_workers > 1 else ([(0, len(labels))] if rows else [])
    pool = ProcessPoolExecutor(n_workers) if n_workers > 1 and len(ranges) > 1 else None
    shared: list[SharedMemory] = []
    try:
        stack: list[np.ndarray | _Shared] = [labels, *mois]
        if pool is not None:  # the stack in shared memory, once, for all workers
            for k, a in enumerate(stack):
                shm, ref = _share(np.ascontiguousarray(a))
                shared.append(shm)
                stack[k] = ref
        labels_in, mois_in = stack[0], stack[1:]

        def run(fn: Any, jobs: list[tuple[Any, ...]]) -> list[Any]:
            nonlocal pool
            if pool is not None:
                try:
                    return list(pool.map(fn, *zip(*jobs, strict=True)))
                except (BrokenProcessPool, OSError) as e:  # e.g. a script without a __main__ guard, or memory
                    log.warning("Worker processes failed (%s): the membranes are measured in one process", e)
                    pool.shutdown(cancel_futures=True)
                    pool = None
            return [fn(*job) for job in jobs]

        jobs1 = []
        for t0, t1 in ranges:
            part = [(t, lab) for t, lab in rows if t0 <= t < t1]
            types = {k: v for k, v in type_of.items() if t0 <= k[0] < t1}
            jobs1.append((labels_in, mois_in, t0, t1, part, types, width, reach, half, list(backgrounds), names,
                          distance))
        first = run(_phase1, jobs1)
        geometry = [c for part in first for c in part[0]]
        places = [p for part in first for p in part[1]]
        frees = [f for part in first for f in part[2]]
        location = locate(geometry, places, frees, names, distance, levels) if geometry else None
        if location is not None:
            location.mode = centre
        lines = []
        for index, c in enumerate(geometry):
            line = c.line
            if location is not None and location.moved:
                shift = location.shift(index, c.typ, c.where)
                if shift.any():  # smoothed along the membrane, so that the band does not jump at a contact's end
                    shift = ndi.uniform_filter1d(shift, size=2 * reach + 1, mode="wrap")
                    line = line + shift[:, None] * c.outward
            lines.append(line)
        jobs2, at = [], 0
        for (t0, t1), part in zip(ranges, first, strict=True):
            n = len(part[0])
            jobs2.append((part[0], lines[at : at + n], mois_in, t0, t1, list(backgrounds), width, inside_only))
            at += n
        second = run(_phase2, jobs2)
    finally:
        if pool is not None:
            pool.shutdown()
        for shm in shared:
            attached = _ATTACHED.pop(shm.name, None)  # attached here too after a fallback to one process
            if attached is not None:
                del attached
            try:
                shm.close()
                shm.unlink()
            except (BufferError, OSError) as e:  # a view still in use: freed when the process ends
                log.debug("Shared memory %s: %s", shm.name, e)
    tracks = [c.track for c in geometry]
    placed = [line + c.corner for c, line in zip(geometry, lines, strict=True)]
    kept = np.array([c.label for c in geometry], dtype=np.int64)
    out = []
    for k in range(len(mois)):
        out.append(Profiles(kept, [p for part in second for p in part[0][k]], tracks,
                            [s for part in second for s in part[1][k]], placed, width, location))
    return out


def selector_problem(sel: str, markers: list[str]) -> str:
    """Why `sel` is not a selection of select() for a panel with these subtype markers ("" when it is one): a cell
    type, 'Immune negative', 'Immune <m>+ <n>+ ...' (subtype markers, in the panel's order, as the tables name the
    subtypes), or '<m>+' (with two or more subtype markers, also with one: every immune cell positive for m)."""
    if sel in ("Tumor", "Immune", "Artefact", "Immune negative"):
        return ""
    if sel.endswith("+") and sel[:-1] in markers and not sel.startswith("Immune "):
        return ""
    if sel.startswith("Immune "):
        rest = sel[len("Immune "):].strip()
        words = _marker_words(rest, markers)  # marker names may contain spaces ("Subtype A+")
        if words is not None and len(set(words)) == len(words):
            ordered = [f"{m}+" for m in markers if f"{m}+" in words]
            if words == ordered:
                return ""
            return f"'{sel}' is not a subset as the tables name it: did you mean 'Immune {' '.join(ordered)}'?"
        if rest.rstrip("+") in markers:
            m = rest.rstrip("+")
            return (f"'{sel}' is not a subset: did you mean '{m}+' (every {m}-positive immune cell) or 'Immune {m}+' "
                    f"({m} positive only)?")
        known = ", ".join(f"'{m}+'" for m in markers) or "none"
        return (f"'{sel}' is not a subset: after 'Immune ' come subtype markers with '+' (of the panel: {known}), "
                "or 'negative'")
    if sel.rstrip("+") in markers:
        m = sel.rstrip("+")
        return f"'{sel}' is not a selection: did you mean '{m}+' (every {m}-positive immune cell)?"
    return f"'{sel}' is not a selection: use Tumor, Immune, 'Immune <subtype>' or '<marker>+'"


def _marker_words(text: str, markers: list[str]) -> list[str] | None:
    """'A+ B c+' -> ['A+', 'B c+'] when every piece is '<marker>+' (longest marker name first), else None."""
    out = []
    names = sorted(markers, key=len, reverse=True)
    while text:
        hit = next((f"{m}+" for m in names if text == f"{m}+" or text.startswith(f"{m}+ ")), None)
        if hit is None:
            return None
        out.append(hit)
        text = text[len(hit):].lstrip()
    return out or None


def selector_problems(settings: MembraneSettings, markers: list[str],
                      markers_of_interest: list[str] | None = None) -> list[str]:
    """The problems of membrane.cells, membrane.partners and membrane.pairs (selector_problem), as messages; with
    markers_of_interest (their names), also the entries of pairs for a name that is not one of them."""
    sides = [("membrane.cells", settings.cells), ("membrane.partners", settings.partners)]
    for marker, entry in settings.pairs.items():
        sides += [(f"membrane.pairs.{marker}.{k}", v) for k, v in entry.items()]
    out = [f"{key}: {msg}" for key, sels in sides for sel in sels if (msg := selector_problem(sel, markers))]
    if markers_of_interest is not None:
        known = ", ".join(f"'{m}'" for m in markers_of_interest) or "none"
        out += [f"membrane.pairs.{m}: '{m}' is not a marker of interest (markers of interest: {known})"
                for m in settings.pairs if m not in markers_of_interest]
    return out


def select(cells: pd.DataFrame, selectors: list[str]) -> pd.Series:
    """Cells matching any of the selectors: a cell type ('Tumor', 'Immune', 'Artefact'), an immune subtype as
    named in the tables ('Immune A+ B+', 'Immune negative'), or a marker ('A+': the immune cells positive for A,
    whatever their other markers)."""
    labels = type_label(cells)
    hit = pd.Series(False, index=cells.index)
    for sel in selectors:
        if sel in ("Tumor", "Immune", "Artefact"):
            hit |= cells["Cell type"] == sel
        elif sel.startswith("Immune "):
            hit |= labels == sel
        elif sel.endswith("+") and sel in cells:
            hit |= (cells["Cell type"] == "Immune") & (cells[sel] == 1)
        else:
            raise ValueError(f"Unknown cell selection '{sel}' (a cell type, 'Immune <subtype>' or '<marker>+')")
    return hit


def per_cell_intensity(cells: pd.DataFrame, membranes: Profiles | None, moi_column: str) -> pd.Series:
    """The marker per cell (by label), background-subtracted: the mean of the membrane that touches no cell, or
    without membrane profiles the mean of the whole cell."""
    if membranes is None:
        return cells.set_index("label")[moi_column]
    return membranes.free()


def carrier(cells: pd.DataFrame, membranes: Profiles | None, moi_column: str) -> tuple[str, dict[str, float]]:
    """The cell type (Tumor or Immune) that carries the marker: the one with the higher median of
    per_cell_intensity (the measure of expression_check), and the medians of both. "" when neither type has cells
    with a value."""
    values = per_cell_intensity(cells, membranes, moi_column)
    medians = {}
    for kind in ("Tumor", "Immune"):
        v = values.reindex(cells.loc[cells["Cell type"] == kind, "label"].astype(int)).to_numpy(dtype=float)
        medians[kind] = float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan")
    finite = {k: v for k, v in medians.items() if np.isfinite(v)}
    return (max(finite, key=lambda k: finite[k]) if finite else ""), medians


def expression_check(cells: pd.DataFrame, membranes: Profiles | None, moi_column: str, wanted: list[str],
                     partners: list[str]) -> dict[str, Any]:
    """Is the marker of interest expressed in the measured cells rather than in their partners? If the partners
    carry it too, their membrane lies against the measured one at the interface and looks like enrichment.

    Compares the median membrane intensity (the membrane that touches no cell; without membrane profiles the mean
    of the whole cell) of the measured cells and of the partners, background-subtracted."""
    values = per_cell_intensity(cells, membranes, moi_column)
    labels = cells["label"].astype(int)

    def median(selectors: list[str]) -> float:
        chosen = values.reindex(labels[select(cells, selectors)]).to_numpy(dtype=float)
        return float(np.nanmedian(chosen)) if np.isfinite(chosen).any() else float("nan")

    a, b = median(wanted), median(partners)
    ratio = b / a if a > 0 else float("inf")
    # the brightest cell type (tumor cells, or one immune subtype): measured cells far below it hardly carry the marker
    by_type = values.reindex(labels).set_axis(type_label(cells).to_numpy()).groupby(level=0).median()
    brightest = float(by_type.max()) if by_type.notna().any() else float("nan")
    if not np.isfinite(a) or not np.isfinite(b):
        verdict = "not checked (no cells of one of the types)"
    elif ratio >= EXPRESSION_SIMILAR[1]:
        verdict = (f"WARNING: the marker of interest is brighter on the partners ({' / '.join(partners)}) than on the "
                   f"measured cells ({' / '.join(wanted)}): check the pair (membrane.cells, membrane.partners; "
                   "empty: the cell type that carries the marker). The interface would show the partners' membrane.")
    elif np.isfinite(brightest) and brightest > 0 and a < EXPRESSION_LOW * brightest:
        verdict = (f"NOT INFORMATIVE: the measured cells ({' / '.join(wanted)}) hardly carry the marker of interest "
                   f"({a / brightest:.0%} of the brightest cell type, {by_type.idxmax()}); their enrichment says "
                   "little.")
    elif ratio >= EXPRESSION_SIMILAR[0]:
        verdict = (f"WARNING: the marker of interest is on both cells at a similar level ({ratio:.0%}): at the "
                   "interface two stained membranes lie side by side, which looks like enrichment. Compare with "
                   "contacts between cells of the same type (membrane.partners).")
    elif ratio >= EXPRESSION_CAUTION:
        verdict = (f"CAUTION: the partners carry {ratio:.0%} of the measured cells' membrane intensity; part of the "
                   "signal at the interface may be theirs (membrane.inside_only limits the band to the measured cell).")
    else:
        verdict = "ok: the marker of interest is mainly on the measured cells"
    level = next((name for prefix, name in CHECK_LEVELS if verdict.startswith(prefix)), "none")
    source = "membrane that touches no cell" if membranes is not None else "whole cell"
    return {"measured cells": " / ".join(wanted), "partners": " / ".join(partners), "measured on": source,
            "level": level,
            "median (measured cells)": a, "median (partners)": b, "partners / measured": ratio, "check": verdict}


# the marker check in a few words (the overview, the pair table), by the start of its verdict
CHECK_LEVELS = [("ok", "mainly on the measured cells"), ("CAUTION", "partners carry part"),
                ("WARNING: the marker of interest is on both", "on both cells"),
                ("WARNING", "brighter on the partners"), ("NOT INFORMATIVE", "measured cells hardly carry it")]
EXPRESSION_CAUTION = 0.25  # partners / measured cells above this: part of the interface signal may be the partners'
EXPRESSION_SIMILAR = (0.8, 1.25)  # partners / measured in this range: on both cells; above: brighter on the partners
EXPRESSION_LOW = 0.25  # measured cells below this fraction of the brightest cell type: they hardly carry the marker


def pair_selections(cells: pd.DataFrame) -> tuple[list[str], list[str]]:
    """The selections of the pair table: rows (membrane of) the cell types, the immune subtypes present and, with two
    or more subtype markers, every marker ('A+'); columns (interface with) the cell types and the subtypes."""
    labels = type_label(cells)
    subtypes = sorted({t for t in labels if t.startswith("Immune ")},
                      key=lambda t: (t != "Immune negative", t.count("+"), t))
    markers = [c for c in cells.columns if c.endswith("+") and not c.startswith(("Immune", "Mean"))]
    rows = ["Tumor", "Immune", *subtypes, *(markers if len(markers) >= 2 else [])]
    return rows, ["Tumor", "Immune", *subtypes]


def pair_table(cells: pd.DataFrame, membranes: Profiles, moi_column: str,
               run_pair: tuple[list[str], list[str]] | None = None) -> pd.DataFrame:
    """Every pair of selections (pair_selections), as the dashboard's 'All pairs': the measured cells with an interface
    with the partners, their median SD above mean and median ratio interface / free membrane (scored as
    analyze_membranes, without the controls), and the marker check of the pair. run_pair (cells, partners): the pair
    of the run, added when it is not one of these (a marker as partner, or several selections), so that the
    comparison of conditions finds it too."""
    rows, cols = pair_selections(cells)
    combos = [([r], [c]) for r in rows for c in cols]
    if run_pair is not None and (list(run_pair[0]), list(run_pair[1])) not in combos:
        combos.append((list(run_pair[0]), list(run_pair[1])))
    labels = cells["label"].astype(int)
    index = {int(lab): k for k, lab in enumerate(membranes.labels)}
    # every position of every membrane at once (2026-11-04: per cell and pair with np.isin took most of a run)
    n = len(membranes.labels)
    lengths = np.array([len(x) for x in membranes.profiles], dtype=np.int64)
    seg = np.repeat(np.arange(n), lengths)
    values = np.concatenate(membranes.profiles).astype(np.float64) if n else np.zeros(0)
    track = np.concatenate(membranes.tracks).astype(np.int64) if n else np.zeros(0, dtype=np.int64)
    free = track == 0
    n_free = np.bincount(seg, free, n)
    sum_free = np.bincount(seg, np.where(free, values, 0.0), n)
    mean_free = np.divide(sum_free, n_free, out=np.full(n, np.nan), where=n_free > 0)
    dev = np.where(free, values - mean_free[seg], 0.0)
    sd_free = np.sqrt(np.divide(np.bincount(seg, dev * dev, n), n_free - 1, out=np.full(n, np.nan),
                                where=n_free > 1))
    top = max(int(track.max(initial=0)), int(labels.to_numpy().max(initial=0))) + 1
    free_of = np.full(top, np.nan)  # every cell's free mean, by label (free_means)
    free_of[membranes.labels.astype(np.int64)] = mean_free
    partner_free = free_of[np.clip(track, 0, None)]
    same_of: dict[str, tuple[float, float]] = {}

    def _same_reference(rs: list[str]) -> tuple[float, float]:
        """The median score and the median log2 ratio of these cells at contacts with cells of the same selection
        (SAME, LOG2_SAME)."""
        own = labels[select(cells, rs)].to_numpy(dtype=np.int64)
        idx = np.array([index[int(lab)] for lab in own if int(lab) in index], dtype=np.int64)
        member = np.zeros(top, dtype=bool)
        member[own] = True
        member[0] = False
        f = member[np.clip(track, 0, None)] & (track > 0)
        n_s = np.bincount(seg, f, n)
        m_s = np.divide(np.bincount(seg, np.where(f, values, 0.0), n), n_s, out=np.full(n, np.nan), where=n_s > 0)
        ok = idx[(n_s[idx] > 0) & (n_free[idx] > 1)]
        sc = np.divide(m_s[ok] - mean_free[ok], sd_free[ok], out=np.full(len(ok), np.nan), where=sd_free[ok] > 0)
        sc = sc[np.isfinite(sc)]
        lr = _log2_of(m_s[ok], mean_free[ok])
        lr = lr[np.isfinite(lr)]
        return (float(np.median(sc)) if len(sc) >= MIN_SAME else float("nan"),
                float(np.median(lr)) if len(lr) >= MIN_SAME else float("nan"))

    out = []
    for rs, cs in combos:
        r, c = " / ".join(rs), " / ".join(cs)
        chosen = np.array([index[int(lab)] for lab in labels[select(cells, rs)] if int(lab) in index], dtype=np.int64)
        is_partner = np.zeros(top, dtype=bool)
        is_partner[labels[select(cells, cs)].to_numpy(dtype=np.int64)] = True
        is_partner[0] = False
        f = is_partner[np.clip(track, 0, None)] & (track > 0)
        n_f = np.bincount(seg, f, n)
        mean_f = np.divide(np.bincount(seg, np.where(f, values, 0.0), n), n_f, out=np.full(n, np.nan),
                           where=n_f > 0)
        known = f & np.isfinite(partner_free)  # partner_excess: the partners' free means over the interface
        n_known = np.bincount(seg, known, n)
        pm = np.divide(np.bincount(seg, np.where(known, partner_free, 0.0), n), n_known, out=np.full(n, np.nan),
                       where=n_known > 0)
        ok = chosen[(n_f[chosen] > 0) & (n_free[chosen] > 1)]
        sd = sd_free[ok]
        scores = np.divide(mean_f[ok] - mean_free[ok], sd, out=np.full(len(ok), np.nan), where=sd > 0)
        same, same_log2 = same_of.setdefault(r, _same_reference(rs))
        log2s = _log2_of(mean_f[ok], mean_free[ok])
        excesses = np.divide(scores * sd - pm[ok], sd, out=np.full(len(ok), np.nan), where=sd > 0)
        ratios = np.divide(mean_f[ok], mean_free[ok], out=np.full(len(ok), np.nan), where=mean_free[ok] != 0)
        check = expression_check(cells, membranes, moi_column, rs, cs)
        out.append({"membrane of": r, "interface with": c, "cells with an interface": len(ok),
                    "median log2 ratio interface/free": _nanmedian(log2s),
                    "median log2 ratio above same-type contacts": _nanmedian(log2s) - same_log2,
                    "median stdDevs above mean": float(np.nanmedian(scores)) if np.isfinite(scores).any() else np.nan,
                    "median ratio interface/free": float(np.nanmedian(ratios)) if np.isfinite(ratios).any() else np.nan,
                    "median excess above both membranes": (float(np.nanmedian(excesses))
                                                           if np.isfinite(excesses).any() else np.nan),
                    "median above same-type contacts": (float(np.nanmedian(scores)) - same
                                                        if np.isfinite(scores).any() else np.nan),
                    "marker check": check["level"], "partners / measured": check["partners / measured"]})
    return pd.DataFrame(out)


def _log2_of(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """log2(a / b) per element; NaN where either is not above 0 (LOG2)."""
    ok = (a > 0) & (b > 0)
    out: np.ndarray = np.full(len(a), np.nan)
    out[ok] = np.log2(a[ok] / b[ok])
    return out


def _nanmedian(v: np.ndarray) -> float:
    return float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan")


def pair_of(settings: MembraneSettings, carried_by: str, marker: str = "") -> tuple[list[str], list[str]]:
    """The cells and partners of the enrichment of one marker, each side on its own: its entry in membrane.pairs,
    else membrane.cells / partners, else automatic: the cell type that carries the marker ("Tumor" or "Immune", see
    carrier; "" counts as Tumor), and the other type."""
    own = "Immune" if carried_by == "Immune" else "Tumor"
    other = "Tumor" if own == "Immune" else "Immune"
    cells, partners = settings.pair_settings(marker)
    return cells or [own], partners or [other]


def analyze_membranes(cells: pd.DataFrame, membranes: Profiles, carried_by: str, settings: MembraneSettings,
                      rng: np.random.Generator | None = None, marker: str = "") -> MembraneResult:
    """The enrichment of the marker of interest at the interface of the chosen cells with their partners
    (pair_of: carried_by is the cell type that carries the marker; marker: its name, for membrane.pairs)."""
    rng = rng or np.random.default_rng()
    wanted, partners = pair_of(settings, carried_by, marker)
    free_of = free_means(membranes)
    chosen = set(cells.loc[select(cells, wanted), "label"].astype(int))
    partner_labels = np.array(sorted(cells.loc[select(cells, partners), "label"].astype(int)), dtype=np.int64)
    frame_of = dict(zip(cells["label"].astype(int), cells["frame"].astype(int), strict=True))
    raw_profiles, iface_flags, free_flags, rows, tracks = [], [], [], [], []
    for label, p, track in zip(membranes.labels, membranes.profiles, membranes.tracks, strict=True):
        if int(label) in chosen:
            raw_profiles.append(p)
            tracks.append(track)
            iface_flags.append(np.isin(track, partner_labels))
            free_flags.append(track == 0)
            rows.append((frame_of[int(label)] + 1, int(label)))
    names = f"{' / '.join(wanted)} cells with {' / '.join(partners)}"
    # the same-type reference: the measured cells' scores at contacts with cells of their own selection
    own_labels = np.array(sorted(chosen), dtype=np.int64)
    same_scores = [s for s in (_score(p[np.isin(t, own_labels)], p[t == 0]) for p, t in zip(raw_profiles, tracks,
                                                                                          strict=True)
                               if np.isin(t, own_labels).any()) if np.isfinite(s)]
    same_ref = float(np.median(same_scores)) if len(same_scores) >= MIN_SAME else float("nan")
    same_log2 = [s for s in (_log2(p[np.isin(t, own_labels)], p[t == 0]) for p, t in zip(raw_profiles, tracks,
                                                                                         strict=True)
                             if np.isin(t, own_labels).any()) if np.isfinite(s)]
    same_log2_ref = float(np.median(same_log2)) if len(same_log2) >= MIN_SAME else float("nan")
    log.info("Same-type contacts of the %s cells: %s", " / ".join(wanted),
             f"median {same_ref:.2f} SD above mean ({len(same_scores)} cells): the reference of '{SAME}'"
             if np.isfinite(same_ref) else f"{len(same_scores)} cells, fewer than {MIN_SAME}: no reference")

    n = len(raw_profiles)
    lengths = np.array([len(p) for p in raw_profiles], dtype=int)
    has = np.array([f.any() for f in iface_flags], dtype=bool)
    iface_lengths = np.array([f.sum() for f in iface_flags if f.any()])
    window = int(round(np.median(iface_lengths))) if len(iface_lengths) else 25
    log.info("Membranes of %d %s, %d with an interface; control window %d px", n, names, has.sum(), window)

    table, aligned, aligned_excess, aligned_log2, control_all, control_log2 = [], [], [], [], [], []
    for p, f, free, track, (frame, label) in zip(raw_profiles, iface_flags, free_flags, tracks, rows, strict=True):
        rest = p[free]  # the membrane that touches no cell
        i_mean = float(p[f].mean()) if f.any() else float("nan")
        r_mean = float(rest.mean()) if len(rest) else float("nan")
        r_sd = float(rest.std(ddof=1)) if len(rest) > 1 else float("nan")
        row = {"Frame": frame, "Label": label, "membrane length": len(p), "interface length": int(f.sum()),
               "interface Mean": i_mean, "non-interface Mean": r_mean, "non-interface stdDev": r_sd,
               "Ratio interface/non-interface": i_mean / r_mean if r_mean else float("nan"),
               LOG2: _log2(p[f], rest), "control log2 ratio": float("nan"),
               "stdDevs above mean": _score(p[f], rest), "control stdDevs above mean": float("nan")}
        partner_mean, sub = partner_excess(p, f, track, free_of)
        ok_sd = bool(r_sd) and np.isfinite(r_sd)
        row[EXCESS] = row["stdDevs above mean"] - partner_mean / r_sd if ok_sd else float("nan")
        row[SAME] = row["stdDevs above mean"] - same_ref
        row[LOG2_SAME] = row[LOG2] - same_log2_ref
        normalized = (p - r_mean) / r_sd if ok_sd else np.full(len(p), np.nan)
        fold = np.full(len(p), np.nan)  # log2 of every position over the free mean (the log2 average profile)
        if r_mean > 0:
            pos = p > 0
            fold[pos] = np.log2(p[pos] / r_mean)
        excess = normalized - sub / r_sd if ok_sd else normalized
        if f.any():
            start, run = max(_runs(f), key=lambda r: r[1])  # the longest interface (the macro: the first one)
            normalized = np.roll(normalized, INTERFACE_OFFSET - (start + run // 2))
            excess = np.roll(excess, INTERFACE_OFFSET - (start + run // 2))
            fold = np.roll(fold, INTERFACE_OFFSET - (start + run // 2))
        elif free.sum() > window:  # windows on the membrane that touches no cell, read as one closed line
            q = p[free]
            scores, folds = [], []
            for s in rng.integers(0, len(q), settings.control_repeats):
                mask = np.zeros(len(q), bool)
                mask[(s + np.arange(window)) % len(q)] = True
                scores.append(_score(q[mask], q[~mask]))
                folds.append(_log2(q[mask], q[~mask]))
            control_all += scores
            control_log2 += folds
            row["control stdDevs above mean"] = _nanmean(scores)
            row["control log2 ratio"] = _nanmean(folds)
        table.append(row)
        aligned.append(normalized)
        aligned_excess.append(excess)
        aligned_log2.append(fold)

    profiles = np.full((n, int(lengths.max(initial=0))), np.nan)
    excess_profiles = profiles.copy()
    log2_profiles = profiles.copy()
    for i, (a, e, g) in enumerate(zip(aligned, aligned_excess, aligned_log2, strict=True)):
        profiles[i, : len(a)] = a
        excess_profiles[i, : len(e)] = e
        log2_profiles[i, : len(g)] = g
    table_df = pd.DataFrame(table)
    with_if = profiles[has]
    length = 0
    if len(with_if) > 1:  # macro: up to round(mean + sd) of the lengths; here of the cells with an interface
        length = min(int(round(lengths[has].mean() + lengths[has].std(ddof=1))), profiles.shape[1])
    def mean_profile(rows_: np.ndarray) -> np.ndarray:
        counts = np.isfinite(rows_[:, :length]).sum(axis=0)
        out: np.ndarray = np.where(counts > 0, np.nansum(rows_[:, :length], axis=0) / np.maximum(counts, 1), np.nan)
        return out

    average = pd.DataFrame({"position (px)": np.arange(length), "profile": mean_profile(with_if),
                            "excess profile": mean_profile(excess_profiles[has]),
                            "log2 ratio profile": mean_profile(log2_profiles[has])})
    control = np.array(control_all, dtype=float)
    sd_values = table_df["stdDevs above mean"].to_numpy(dtype=float) if n else np.zeros(0)
    top, _, n_bins = histogram_scale(sd_values, control, settings.max_histogram_value, settings.histogram_bin_width)
    return MembraneResult(table_df, profiles, lengths, has, average, _histogram(sd_values, top, n_bins),
                          _histogram(control, top, n_bins), control, window, excess_profiles,
                          np.array(control_log2, dtype=float))


def _nanmean(v: list[float]) -> float:
    a = np.asarray(v, dtype=float)
    return float(np.nanmean(a)) if np.isfinite(a).any() else float("nan")
