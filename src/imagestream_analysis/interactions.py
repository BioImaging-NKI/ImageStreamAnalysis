"""Which cells touch which, within each event (ClassifyCellTypes 2.8, lines 479-484 and 609-632; replaces the
Summarize_Interactions notebook).

2.8 uses MorphoLibJ's Region Adjacency Graph on the eroded labels. Here two cells interact when one comes within
`reach` pixels (4-connected steps) of the other: the same reach as the interfaces of the membrane analysis, so that
a cell with an interface always has the partner in its list. Cells are named by their label, which is unique in
the dataset (the 'Label' column of the cell table)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def touching_pairs(labels: np.ndarray, reach: int) -> np.ndarray:
    """(P, 2) label pairs (a < b) of cells within `reach` pixels of each other, in the same event."""
    t, h, w = labels.shape
    pairs = set()
    for i, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        ys = slice(max(sl[1].start - reach, 0), min(sl[1].stop + reach, h))
        xs = slice(max(sl[2].start - reach, 0), min(sl[2].stop + reach, w))
        sub = labels[sl[0].start, ys, xs]
        near = ndi.distance_transform_cdt(sub != i, metric="taxicab") <= reach
        for j in np.unique(sub[near & (sub != i) & (sub > 0)]):
            pairs.add((min(i, int(j)), max(i, int(j))))
    return np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)


def type_label(cells: pd.DataFrame) -> pd.Series:
    """The cell type as one name: 'Tumor', 'Artefact', 'Immune', or 'Immune <subtype>'."""
    subtype = cells["Subtype"].fillna("") if "Subtype" in cells else pd.Series("", index=cells.index)
    return cells["Cell type"].where(subtype == "", cells["Cell type"] + " " + subtype)


def add_neighbours(cells: pd.DataFrame, pairs: np.ndarray) -> pd.DataFrame:
    """The cell table with 'Neighbours' (number), 'Touching (labels)' and 'Neighbour types'."""
    out = cells.copy()
    partners: dict[int, list[int]] = {}
    for a, b in pairs:
        partners.setdefault(int(a), []).append(int(b))
        partners.setdefault(int(b), []).append(int(a))
    type_of = dict(zip(cells["label"].astype(int), type_label(cells), strict=True))
    lists = [sorted(partners.get(int(lab), [])) for lab in cells["label"]]
    out["Neighbours"] = [len(p) for p in lists]
    out["Touching (labels)"] = [",".join(map(str, p)) for p in lists]
    out["Neighbour types"] = ["; ".join(sorted(type_of[j] for j in p)) for p in lists]
    return out


def contact_table(cells: pd.DataFrame, pairs: np.ndarray) -> pd.DataFrame:
    """Per cell type and neighbour type: the contacts (pairs), the cells of the type, and how many of them touch at
    least one cell of the neighbour type. Every pair counts from both sides (A–B and B–A)."""
    types = type_label(cells)
    type_of = dict(zip(cells["label"].astype(int), types, strict=True))
    rows = []
    for a, b in pairs:
        ta, tb = type_of[int(a)], type_of[int(b)]
        rows += [(ta, tb, int(a)), (tb, ta, int(b))]
    contacts = pd.DataFrame(rows, columns=["cell type", "neighbour type", "label"])
    n_cells = types.value_counts()
    table = (contacts.groupby(["cell type", "neighbour type"])
             .agg(contacts=("label", "size"), **{"cells touching": ("label", "nunique")}).reset_index())
    table.insert(2, "cells", table["cell type"].map(n_cells).astype(int))
    table["fraction touching"] = table["cells touching"] / table["cells"]
    return table.sort_values(["cell type", "neighbour type"]).reset_index(drop=True)


def multiplets(cells: pd.DataFrame, pairs: np.ndarray, n_events: int) -> pd.DataFrame:
    """The number of events with n cells and with n interactions (2.8: the 'multiplets' table)."""
    per_event_cells = np.bincount(cells["frame"].to_numpy(dtype=int), minlength=n_events)
    per_event_pairs = interactions_per_event(cells, pairs, n_events)
    top = int(max(per_event_cells.max(initial=0), per_event_pairs.max(initial=0)))
    return pd.DataFrame({"n": np.arange(top + 1),
                         "events with n cells": np.bincount(per_event_cells, minlength=top + 1),
                         "events with n interactions": np.bincount(per_event_pairs, minlength=top + 1)})


def interactions_per_event(cells: pd.DataFrame, pairs: np.ndarray, n_events: int) -> np.ndarray:
    """The number of interacting pairs in each event."""
    frame_of = dict(zip(cells["label"].astype(int), cells["frame"].astype(int), strict=True))
    frames = np.array([frame_of[int(a)] for a, _ in pairs], dtype=np.int64)
    out: np.ndarray = np.bincount(frames, minlength=n_events).astype(int)
    return out
