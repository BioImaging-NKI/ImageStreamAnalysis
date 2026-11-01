"""Membranes, interfaces and the per-event measurements of the marker of interest (macro lines 509-693)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import imagej
from .classification import IMMUNE, TUMOR


@dataclass
class InterfaceMasks:
    """Binary (T, Y, X) masks."""

    immune_membrane: np.ndarray
    tumor_membrane: np.ndarray
    immune_interface: np.ndarray
    tumor_interface: np.ndarray
    immune_not_in_interface: np.ndarray
    tumor_not_in_interface: np.ndarray


def interface_masks(classes: np.ndarray, thickness: int) -> InterfaceMasks:
    """Membrane = the outer `thickness` pixels of each class mask; interface of an immune cell = its pixels within
    `thickness` pixels of a tumor cell (and vice versa); 'not in interface' = interface XOR membrane, as the macro;
    since the interface lies inside the membrane band this is membrane AND NOT interface (SPEC, I6)."""
    t_mask = classes == IMMUNE
    tu_mask = classes == TUMOR
    t_membrane = t_mask & ~imagej.erode(t_mask, thickness)
    tu_membrane = tu_mask & ~imagej.erode(tu_mask, thickness)
    t_interface = t_mask & imagej.dilate(tu_mask, thickness)
    tu_interface = tu_mask & imagej.dilate(t_mask, thickness)
    return InterfaceMasks(t_membrane, tu_membrane, t_interface, tu_interface,
                          t_interface ^ t_membrane, tu_interface ^ tu_membrane)


def _masked_means(image: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Mean of image within mask, per event; NaN for an empty mask (CLIJ2 getMeanOfMaskedPixels)."""
    t = image.shape[0]
    m = mask.reshape(t, -1)
    total = np.where(m, image.reshape(t, -1), 0).sum(axis=1, dtype=np.float64)
    count = m.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        means: np.ndarray = total / count
    return means


MOI_COLUMNS = ["mean immune interface", "mean immune NOT in interface", "mean tumor interface",
               "mean tumor NOT in interface", "Immune membrane ratio", "Tumor membrane ratio"]


def event_table(cells: pd.DataFrame, n_events: int, masks: InterfaceMasks, moi: np.ndarray | None) -> pd.DataFrame:
    """The 'All cell combinations' table: one row per event (Frame is 1-based). Without a marker of interest
    (moi None) the marker columns (MOI_COLUMNS) are left out.

    The macro writes a row only for events with at least one cell, but fills the measurements by event index,
    so its rows shift when an event has no cells (issue I5); here every event has its own row."""
    counts = pd.crosstab(cells["frame"], cells["Cell type"]).reindex(range(n_events), fill_value=0)
    values = moi if moi is not None else np.zeros(masks.immune_interface.shape, dtype=np.float32)
    df = pd.DataFrame({
        "Frame": np.arange(1, n_events + 1),
        "Immune cells": counts.get("Immune", pd.Series(0, index=counts.index)).to_numpy(),
        "Tumor cells": counts.get("Tumor", pd.Series(0, index=counts.index)).to_numpy(),
        "mean immune interface": _masked_means(values, masks.immune_interface),
        "mean immune NOT in interface": _masked_means(values, masks.immune_not_in_interface),
        "mean tumor interface": _masked_means(values, masks.tumor_interface),
        "mean tumor NOT in interface": _masked_means(values, masks.tumor_not_in_interface),
    })
    with np.errstate(divide="ignore", invalid="ignore"):
        df["Immune membrane ratio"] = df["mean immune interface"] / df["mean immune NOT in interface"]
        df["Tumor membrane ratio"] = df["mean tumor interface"] / df["mean tumor NOT in interface"]
    # whether a tumor and an immune cell touch (the interface mask is not empty), with or without a marker
    df["Tumor-immune interface"] = masks.immune_interface.reshape(n_events, -1).any(axis=1)
    return df if moi is not None else df.drop(columns=MOI_COLUMNS)


def doublets(events: pd.DataFrame) -> pd.DataFrame:
    """Events with exactly one immune cell and one tumor cell that touch ('Single Tumor Single Immune cell')."""
    keep = (events["Immune cells"] == 1) & (events["Tumor cells"] == 1) & events["Tumor-immune interface"]
    return events[keep].reset_index(drop=True)
