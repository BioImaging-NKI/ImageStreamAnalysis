"""Per-cell measurements: the first stage of the analysis, before the classification.

The normalized tumor and T-cell marker intensities of the macro v1.8 (lines 325-370), plus (new, as
ClassifyCellTypes 2.8) the area in µm² and the background-subtracted mean of every loaded channel, and whether a cell
is cut by the edge of its event image."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .io import EventStack

log = logging.getLogger(__name__)


def frames_of(labels: np.ndarray) -> np.ndarray:
    """The event (0-based) of each label 1..N: the first event in which it occurs (-1 if absent)."""
    n = int(labels.max(initial=0))
    frame_of = np.full(n + 1, -1, dtype=np.int64)
    present = labels.reshape(labels.shape[0], -1)
    for t in range(labels.shape[0]):
        vals = np.unique(present[t])
        vals = vals[(vals > 0) & (frame_of[vals] < 0)]
        frame_of[vals] = t
    out: np.ndarray = frame_of[1:]
    return out


def edge_distance(labels: np.ndarray, extents: np.ndarray) -> np.ndarray:
    """For each label 1..N, the distance in pixels from the cell to the edge of its event's own image (0 = on the
    outermost row or column). extents: (T, 4) y0, x0, y1, x1 (EventStack.image_extents)."""
    t, h, w = labels.shape
    yy = np.arange(h, dtype=np.int16)[None, :, None]
    xx = np.arange(w, dtype=np.int16)[None, None, :]
    e = extents.astype(np.int16)
    dist = np.minimum(np.minimum(yy - e[:, 0, None, None], e[:, 2, None, None] - 1 - yy),
                      np.minimum(xx - e[:, 1, None, None], e[:, 3, None, None] - 1 - xx))
    index = np.arange(1, int(labels.max(initial=0)) + 1)
    out: np.ndarray = np.asarray(ndi.minimum(dist, labels, index), dtype=np.int64)
    return out


def edge_labels(labels: np.ndarray, extents: np.ndarray, erosion: int) -> np.ndarray:
    """Labels 1..N -> True for a cell cut by the edge of its event image: within `erosion` pixels of it (the cells
    were shrunk from the background by that much, so a cell that touched the edge now lies `erosion` pixels away)."""
    out: np.ndarray = edge_distance(labels, extents) <= erosion
    return out


def channel_column(number: int, name: str) -> str:
    return f"Mean ch{number} ({name})"


def measure_cells(labels: np.ndarray, stack: EventStack, tumor: np.ndarray, t_cell: np.ndarray,
                  edge: np.ndarray) -> pd.DataFrame:
    """One row per cell, ordered by event and then label (the order of the macro's ROI Manager).

    labels: (T, Y, X) cleaned labels 1..N; tumor, t_cell: the normalized marker channels (T, Y, X); edge: per
    label 1..N, whether the cell is cut by the edge of its event image.

    Columns: label, frame (0-based), T-cell and tumor intensity (mean of the normalized markers), Area (pixels),
    Marker ratio, Marker ratio over Area, Area (µm²), Edge, and 'Mean ch<N> (<name>)' for every loaded channel: the
    mean of the stored values minus the channel's background (median)."""
    n = int(labels.max(initial=0))
    index = np.arange(1, n + 1)
    area = ndi.sum_labels(np.ones_like(labels, dtype=np.int64), labels, index).astype(np.int64)
    df = pd.DataFrame({"label": index, "frame": frames_of(labels),
                       "T-cell intensity": ndi.mean(t_cell, labels, index),
                       "Tumor intensity": ndi.mean(tumor, labels, index), "Area": area})
    with np.errstate(divide="ignore", invalid="ignore"):
        df["Marker ratio"] = df["T-cell intensity"] / df["Tumor intensity"]
        df["Marker ratio over Area"] = df["Marker ratio"] / df["Area"]
    df["Area (µm²)"] = area * stack.pixel_size**2
    df["Edge"] = edge.astype(bool)
    numbers = stack.channel_numbers or list(range(1, stack.data.shape[1] + 1))
    for number, name in zip(numbers, stack.channel_names, strict=True):
        means = ndi.mean(stack.channel(number), labels, index)
        df[channel_column(number, name)] = np.asarray(means, dtype=np.float64) - stack.median(number)
    return df.sort_values(["frame", "label"], kind="stable").reset_index(drop=True)
