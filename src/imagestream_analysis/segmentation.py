"""Cell segmentation: the Cellpose input (macro lines 199-279), Cellpose in a subprocess, and the clean-up of
the labels (lines 281-317)."""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import imagej
from .config import MarkerLocation, Settings
from .io import EventStack, grid_columns, montage_to_stack, read_labelmap, stack_to_montage

log = logging.getLogger(__name__)
RUNNER = Path(__file__).with_name("_cellpose_runner.py")


def normalize_marker(channel: np.ndarray) -> tuple[np.ndarray, float]:
    """The macro's normalizeImage(): divide by the median over the events of the per-event Otsu threshold."""
    thresholds = [imagej.otsu_threshold_dark(frame) for frame in channel]
    factor = imagej.median_upper(thresholds)
    return channel / np.float32(factor), factor


@dataclass
class Markers:
    """The normalized tumor and T-cell marker channels, (T, Y, X) float32, used for segmentation and
    classification."""

    tumor: np.ndarray
    t_cell: np.ndarray
    tumor_factor: float
    t_cell_factor: float


def markers(stack: EventStack, settings: Settings) -> Markers:
    ch = settings.channels
    tumor, f_tumor = normalize_marker(stack.channel(ch.tumor))
    t_cell, f_t = normalize_marker(stack.channel(ch.t_cell))
    log.info("Marker normalization factors (median Otsu threshold): tumor %.4g, T cell %.4g", f_tumor, f_t)
    return Markers(tumor, t_cell, f_tumor, f_t)


def montage_columns(stack: EventStack, settings: Settings) -> int:
    t, _, h, w = stack.data.shape
    columns = settings.segmentation.montage_columns
    return min(columns, t) if columns > 0 else grid_columns(t, h, w)


def cellpose_input(stack: EventStack, m: Markers, settings: Settings) -> tuple[np.ndarray, list[int]]:
    """The image for Cellpose: a montage of all events (a grid, row by row; the macro uses one row: SPEC I2), and
    the Cellpose 3 channels. Empty tiles at the end of the grid are 0.

    nuclear markers: channel 1 = variance (r=2) of the brightfield, channel 2 = tumor + T-cell marker (normalized);
    Cellpose with ch1=1, ch2=2. membrane markers: one channel, tumor + T-cell marker + brightfield variance, each
    normalized; Cellpose with ch1=1, ch2=0 (grayscale)."""
    columns = montage_columns(stack, settings)
    bf_var = imagej.variance_filter(stack.channel(settings.channels.brightfield), 2)
    nuclei = m.tumor + m.t_cell
    if settings.channels.marker_location == MarkerLocation.NUCLEAR:
        image = stack_to_montage(np.stack([bf_var, nuclei], axis=-1), columns)
        return image.astype(np.float32), [1, 2]
    bf_norm, _ = normalize_marker(bf_var)
    return stack_to_montage(nuclei + bf_norm, columns).astype(np.float32), [0, 0]


def run_cellpose(image: np.ndarray, channels: list[int], settings: Settings) -> tuple[np.ndarray, dict[str, Any]]:
    seg = settings.segmentation
    python = str(seg.cellpose_python) if seg.cellpose_python else sys.executable
    params = {"model": seg.model, "diameter": seg.diameter, "flow_threshold": seg.flow_threshold,
              "cellprob_threshold": seg.cellprob_threshold, "use_gpu": seg.use_gpu, "channels": channels}
    with tempfile.TemporaryDirectory(prefix="isa_cellpose_") as tmp:
        src, dst = Path(tmp, "in.npy"), Path(tmp, "out.npy")
        np.save(src, image)
        log.info("Cellpose (%s, model %s, diameter %s) on a %d x %d px image", python, seg.model, seg.diameter,
                 image.shape[1], image.shape[0])
        proc = subprocess.run([python, "-I", str(RUNNER), str(src), str(dst), json.dumps(params)],
                              check=False, capture_output=True, text=True,
                              creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        if proc.returncode != 0:
            raise RuntimeError(f"Cellpose failed (exit code {proc.returncode}):\n{proc.stderr[-3000:]}")
        info = {}
        for line in proc.stdout.splitlines():
            if line.startswith("ISA_INFO "):
                info = json.loads(line[len("ISA_INFO "):])
        log.info("Cellpose done: %s", info)
        return np.load(dst), info


def segment(stack: EventStack, m: Markers, settings: Settings) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    """Raw labels (T, Y, X) as Cellpose made them, plus the montage to save as labelmap_cells (None when loaded
    from disk)."""
    seg = settings.segmentation
    t, _, h, w = stack.data.shape
    if seg.method == "labelmap":
        if not seg.labelmap_file:
            raise ValueError("segmentation.method = 'labelmap' needs segmentation.labelmap_file")
        log.info("Loading labelmap %s", seg.labelmap_file)
        return read_labelmap(seg.labelmap_file, t, h, w), None, {}
    if seg.method != "cellpose":
        raise ValueError(f"Unknown segmentation.method '{seg.method}' (cellpose or labelmap)")
    image, channels = cellpose_input(stack, m, settings)
    montage, info = run_cellpose(image, channels, settings)
    return montage_to_stack(montage, h, w)[:t], montage, info


def split_across_events(labels: np.ndarray) -> np.ndarray:
    """Give a label that occurs in more than one event (a cell that Cellpose drew across the border of two tiles of
    the montage) a separate label in each event (SPEC I15). Labels stay ordered by event, then label."""
    t = labels.shape[0]
    pairs = np.zeros(labels.shape, dtype=np.int64)
    fg = labels > 0
    pairs[fg] = (np.arange(t, dtype=np.int64)[:, None, None] * (int(labels.max(initial=0)) + 1)
                 + labels)[fg]
    values, inverse = np.unique(pairs, return_inverse=True)
    n_split = (len(values) - 1) - len(np.unique(labels[fg]))
    if n_split:
        log.info("%d cell(s) crossed the border between two events in the montage: split per event", n_split)
    out: np.ndarray = inverse.reshape(labels.shape).astype(np.int32)  # 0 stays 0: it is the smallest value
    return out


def clean_labels(labels: np.ndarray, erosion: int, min_size: int, max_size: int) -> np.ndarray:
    """Split cells across events (SPEC I15); close index gaps; shrink the cells from the background by `erosion`
    pixels; remove cells outside the size range; close gaps again. As in the macro, the size filter is only applied
    when erosion > 0 (issue I3)."""
    labels = imagej.close_index_gaps(split_across_events(labels))
    if erosion <= 0:
        return labels
    eroded = labels * imagej.erode(labels > 0, erosion)
    sizes = np.bincount(eroded.ravel())
    keep = (sizes >= min_size) & (sizes <= max_size)
    keep[0] = False
    eroded = np.where(keep[eroded], eroded, 0)
    return imagej.close_index_gaps(eroded)


def remove_labels(labels: np.ndarray, remove: np.ndarray) -> np.ndarray:
    """Remove the labels for which remove[label - 1] is True and close the index gaps."""
    drop = np.concatenate([[False], remove.astype(bool)])
    return imagej.close_index_gaps(np.where(drop[labels], 0, labels))
