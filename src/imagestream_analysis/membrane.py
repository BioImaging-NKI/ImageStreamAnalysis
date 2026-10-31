"""Membrane profiles of the cells of interest and the enrichment at the interface (macro lines 700-996).

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

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure

from .config import CellOfInterest, MembraneSettings

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


def straighten(image: np.ndarray, line: np.ndarray, width: int) -> np.ndarray:
    """Fiji Straighten along a closed line: `width` samples across the line (bilinear), (width, n)."""
    tangent = (np.roll(line, -1, axis=0) - np.roll(line, 1, axis=0)) / 2
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-9)
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    offsets = np.arange(width) - (width - 1) / 2
    pts = line[None] + offsets[:, None, None] * normal[None]
    return ndi.map_coordinates(image, [pts[..., 0], pts[..., 1]], order=1, mode="nearest")


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


def analyze_membranes(labels: np.ndarray, cells: pd.DataFrame, interface: np.ndarray, moi: np.ndarray,
                      background: float, cell_of_interest: CellOfInterest, settings: MembraneSettings,
                      rng: np.random.Generator | None = None) -> MembraneResult:
    """labels (T, Y, X); cells: the classification table; interface: (T, Y, X) interface mask of the class of
    interest; moi: (T, Y, X) raw marker of interest, `background` is subtracted."""
    rng = rng or np.random.default_rng()
    wanted = "Immune" if cell_of_interest == CellOfInterest.IMMUNE_CELL else "Tumor"
    roi = cells[cells["Cell type"] == wanted]
    width = settings.straightening_width
    objects = {t: ndi.find_objects(labels[t]) for t in np.unique(roi["frame"])}
    raw_profiles, iface_flags, rows = [], [], []
    for _, cell in roi.iterrows():
        t, label = int(cell["frame"]), int(cell["label"])
        sl = objects[t][label - 1] if label - 1 < len(objects[t]) else None
        if sl is None:
            continue
        pad = width  # the straightened band reaches outside the cell
        y0, x0 = max(sl[0].start - pad, 0), max(sl[1].start - pad, 0)
        box = (slice(y0, sl[0].stop + pad), slice(x0, sl[1].stop + pad))
        line = outline(labels[t][box] == label)
        if line is None:
            continue
        band = straighten(moi[t][box].astype(np.float64) - background, line, width)
        own_interface = interface[t][box] & (labels[t][box] == label)
        iband = straighten(own_interface.astype(np.float32), line, width)
        raw_profiles.append(band.max(axis=0))  # macro: maximum projection across the membrane
        iface_flags.append(iband.max(axis=0) > 0)
        rows.append((t + 1, label))

    n = len(raw_profiles)
    lengths = np.array([len(p) for p in raw_profiles], dtype=int)
    has = np.array([f.any() for f in iface_flags], dtype=bool)
    iface_lengths = np.array([f.sum() for f in iface_flags if f.any()])
    window = int(round(np.median(iface_lengths))) if len(iface_lengths) else 25
    log.info("Membranes of %d %s cells, %d with an interface; control window %d px", n, wanted, has.sum(), window)

    table, aligned, control_all = [], [], []
    for p, f, (frame, label) in zip(raw_profiles, iface_flags, rows, strict=True):
        rest = p[~f]
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
        elif len(p) > window:
            scores = []
            for s in rng.integers(0, len(p), settings.control_repeats):
                idx = (s + np.arange(window)) % len(p)
                mask = np.zeros(len(p), bool)
                mask[idx] = True
                scores.append(_score(p[mask], p[~mask]))
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
