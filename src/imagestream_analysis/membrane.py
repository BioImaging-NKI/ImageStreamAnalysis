"""Membrane profiles and the enrichment at the interface (macro lines 700-996), between two cell types of choice.

Every cell (artefacts excepted) gets its membrane profile of the marker of interest and a neighbour track: for each
position along the outline, the label of the nearest other cell within the interface reach (0 = none). A pair of
cell types (`membrane.cells`, `membrane.partners`) then only chooses which positions are the interface: those that
touch a partner; the rest of the membrane is the part that touches no cell. (Decided 2026-10-31: with several cell
types the interfaces between all of them are kept by position, and any pair is a selection.)

Deviations from the macro, decided 2026-10-22 (SPEC I8-I11):
- one membrane per cell (the macro: one per blob of touching cells of the same class);
- the length of each profile is kept explicitly (the macro marks the end with NaN, which also hits true zeros);
- the control for cells without an interface: `control_repeats` random windows per cell (unseeded), as long as the
  median interface length (the macro: one window of 25 px);
- the average profile only uses cells with an interface.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure

from .config import CellOfInterest, MembraneSettings
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


def _score(window: np.ndarray, rest: np.ndarray) -> float:
    """(window mean - rest mean) / rest sd (sample sd, as ImageJ)."""
    if len(window) == 0 or len(rest) < 2:
        return float("nan")
    sd = rest.std(ddof=1)
    return float((window.mean() - rest.mean()) / sd) if sd > 0 else float("nan")


def _histogram(values: np.ndarray, settings: MembraneSettings) -> pd.DataFrame:
    """As the macro's Distribution: round((max+5)/binWidth) bins over [-5, max]; values outside are not counted."""
    hi = settings.max_histogram_value
    n_bins = int(round((hi - HISTOGRAM_MIN) / settings.histogram_bin_width))
    counts, edges = np.histogram(values[np.isfinite(values)], bins=n_bins, range=(HISTOGRAM_MIN, hi))
    return pd.DataFrame({"bin_start": edges[:-1], "counts": counts})


@dataclass
class Profiles:
    """The membrane of every cell: its profile (marker of interest minus background, maximum across the membrane)
    and its neighbour track (label of the nearest other cell within the reach, 0 = none), same length."""

    labels: np.ndarray  # the cell (label) of each profile
    profiles: list[np.ndarray]
    tracks: list[np.ndarray]


def membrane_profiles(labels: np.ndarray, cells: pd.DataFrame, moi: np.ndarray, background: float, width: int,
                      reach: int, inside_only: bool = False) -> Profiles:
    """Profiles and neighbour tracks of all cells but the artefacts. labels (T, Y, X); moi: (T, Y, X) raw marker of
    interest, `background` is subtracted; reach: the interface reach (membrane.thickness), in 4-connected steps.

    The straightened band is `width` samples across the outline, half of them outside the cell: at an interface
    they lie in the partner cell. inside_only: the maximum is taken over the samples inside the cell only, so that
    the membrane of a partner that also carries the marker is not measured as this cell's."""
    roi = cells[cells["Cell type"] != "Artefact"]
    objects = {t: ndi.find_objects(labels[t]) for t in np.unique(roi["frame"])}
    out = Profiles(np.zeros(0, dtype=np.int64), [], [])
    kept = []
    for t, label in zip(roi["frame"].astype(int), roi["label"].astype(int), strict=True):
        sl = objects[t][label - 1] if label - 1 < len(objects[t]) else None
        if sl is None:
            continue
        pad = width + reach  # the straightened band reaches outside the cell; neighbours up to `reach` further
        y0, x0 = max(sl[0].start - pad, 0), max(sl[1].start - pad, 0)
        box = (slice(y0, sl[0].stop + pad), slice(x0, sl[1].stop + pad))
        sub = labels[t][box]
        own = sub == label
        line = outline(own)
        if line is None:
            continue
        # every own pixel within `reach` of another cell: the label of the nearest one (as imagej.dilate: taxicab)
        others = np.where(own, 0, sub)
        dist, (iy, ix) = ndi.distance_transform_cdt(others == 0, metric="taxicab", return_indices=True)
        nearest = np.where(own & (dist <= reach), others[iy, ix], 0)
        band = straighten(moi[t][box].astype(np.float64) - background, line, width)
        if inside_only:
            inside = straighten(own.astype(np.float64), line, width, order=0) > 0
            inside[(width - 1) // 2] |= ~inside.any(axis=0)  # a sliver of a cell: at least the middle sample
            band = np.where(inside, band, -np.inf)
        out.profiles.append(band.max(axis=0))  # macro: maximum projection across the membrane
        out.tracks.append(straighten(nearest.astype(np.float64), line, width, order=0).max(axis=0).astype(np.int64))
        kept.append(label)
    out.labels = np.array(kept, dtype=np.int64)
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


def expression_check(cells: pd.DataFrame, membranes: Profiles | None, moi_column: str, wanted: list[str],
                     partners: list[str]) -> dict[str, Any]:
    """Is the marker of interest expressed in the measured cells rather than in their partners? If the partners
    carry it too, their membrane lies against the measured one at the interface and looks like enrichment.

    Compares the median membrane intensity (the membrane that touches no cell; without membrane profiles the mean
    of the whole cell) of the measured cells and of the partners, background-subtracted."""
    def per_cell() -> pd.Series:
        if membranes is None:
            return cells.set_index("label")[moi_column]
        free = [float(p[t == 0].mean()) if (t == 0).any() else np.nan
                for p, t in zip(membranes.profiles, membranes.tracks, strict=True)]
        return pd.Series(free, index=membranes.labels)

    values = per_cell()
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
                   f"measured cells ({' / '.join(wanted)}): check where it is expressed (channels.cell_of_interest, "
                   "membrane.cells). The interface would show the partners' membrane.")
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
    source = "membrane that touches no cell" if membranes is not None else "whole cell"
    return {"measured cells": " / ".join(wanted), "partners": " / ".join(partners), "measured on": source,
            "median (measured cells)": a, "median (partners)": b, "partners / measured": ratio, "check": verdict}


EXPRESSION_CAUTION = 0.25  # partners / measured cells above this: part of the interface signal may be the partners'
EXPRESSION_SIMILAR = (0.8, 1.25)  # partners / measured in this range: on both cells; above: brighter on the partners
EXPRESSION_LOW = 0.25  # measured cells below this fraction of the brightest cell type: they hardly carry the marker


def pair_of(settings: MembraneSettings, cell_of_interest: CellOfInterest) -> tuple[list[str], list[str]]:
    """The cells and partners of the enrichment; empty settings: the cell of interest against the other type."""
    own = "Immune" if cell_of_interest == CellOfInterest.IMMUNE_CELL else "Tumor"
    other = "Tumor" if own == "Immune" else "Immune"
    return list(settings.cells) or [own], list(settings.partners) or [other]


def analyze_membranes(cells: pd.DataFrame, membranes: Profiles, cell_of_interest: CellOfInterest,
                      settings: MembraneSettings, rng: np.random.Generator | None = None) -> MembraneResult:
    """The enrichment of the marker of interest at the interface of the chosen cells with their partners."""
    rng = rng or np.random.default_rng()
    wanted, partners = pair_of(settings, cell_of_interest)
    chosen = set(cells.loc[select(cells, wanted), "label"].astype(int))
    partner_labels = np.array(sorted(cells.loc[select(cells, partners), "label"].astype(int)), dtype=np.int64)
    frame_of = dict(zip(cells["label"].astype(int), cells["frame"].astype(int), strict=True))
    raw_profiles, iface_flags, free_flags, rows = [], [], [], []
    for label, p, track in zip(membranes.labels, membranes.profiles, membranes.tracks, strict=True):
        if int(label) in chosen:
            raw_profiles.append(p)
            iface_flags.append(np.isin(track, partner_labels))
            free_flags.append(track == 0)
            rows.append((frame_of[int(label)] + 1, int(label)))
    names = f"{' / '.join(wanted)} cells with {' / '.join(partners)}"

    n = len(raw_profiles)
    lengths = np.array([len(p) for p in raw_profiles], dtype=int)
    has = np.array([f.any() for f in iface_flags], dtype=bool)
    iface_lengths = np.array([f.sum() for f in iface_flags if f.any()])
    window = int(round(np.median(iface_lengths))) if len(iface_lengths) else 25
    log.info("Membranes of %d %s, %d with an interface; control window %d px", n, names, has.sum(), window)

    table, aligned, control_all = [], [], []
    for p, f, free, (frame, label) in zip(raw_profiles, iface_flags, free_flags, rows, strict=True):
        rest = p[free]  # the membrane that touches no cell
        i_mean = float(p[f].mean()) if f.any() else float("nan")
        r_mean = float(rest.mean()) if len(rest) else float("nan")
        r_sd = float(rest.std(ddof=1)) if len(rest) > 1 else float("nan")
        row = {"Frame": frame, "Label": label, "membrane length": len(p), "interface length": int(f.sum()),
               "interface Mean": i_mean, "non-interface Mean": r_mean, "non-interface stdDev": r_sd,
               "Ratio interface/non-interface": i_mean / r_mean if r_mean else float("nan"),
               "stdDevs above mean": _score(p[f], rest), "control stdDevs above mean": float("nan")}
        normalized = (p - r_mean) / r_sd if r_sd and np.isfinite(r_sd) else np.full(len(p), np.nan)
        if f.any():
            start, run = max(_runs(f), key=lambda r: r[1])  # the longest interface (the macro: the first one)
            normalized = np.roll(normalized, INTERFACE_OFFSET - (start + run // 2))
        elif free.sum() > window:  # windows on the membrane that touches no cell, read as one closed line
            q = p[free]
            scores = []
            for s in rng.integers(0, len(q), settings.control_repeats):
                mask = np.zeros(len(q), bool)
                mask[(s + np.arange(window)) % len(q)] = True
                scores.append(_score(q[mask], q[~mask]))
            control_all += scores
            row["control stdDevs above mean"] = float(np.nanmean(scores))
        table.append(row)
        aligned.append(normalized)

    profiles = np.full((n, int(lengths.max(initial=0))), np.nan)
    for i, a in enumerate(aligned):
        profiles[i, : len(a)] = a
    table_df = pd.DataFrame(table)
    with_if = profiles[has]
    length = 0
    if len(with_if) > 1:  # macro: up to round(mean + sd) of the lengths; here of the cells with an interface
        length = min(int(round(lengths[has].mean() + lengths[has].std(ddof=1))), profiles.shape[1])
    counts = np.isfinite(with_if[:, :length]).sum(axis=0)
    avg = np.where(counts > 0, np.nansum(with_if[:, :length], axis=0) / np.maximum(counts, 1), np.nan)
    average = pd.DataFrame({"position (px)": np.arange(length), "profile": avg})
    control = np.array(control_all, dtype=float)
    return MembraneResult(table_df, profiles, lengths, has, average,
                          _histogram(table_df["stdDevs above mean"].to_numpy(dtype=float) if n else np.zeros(0),
                                     settings),
                          _histogram(control, settings), control, window)
