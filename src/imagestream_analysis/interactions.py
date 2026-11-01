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


MIN_CELLS_FOR_CHANCE = 3  # events with fewer cells cannot show a partner preference (see contact_table)


def contact_table(cells: pd.DataFrame, pairs: np.ndarray) -> pd.DataFrame:
    """Per cell type and neighbour type: the cells of the type, the contacts (every pair counts from both sides:
    A–B and B–A), and how many cells of the type touch at least one of that neighbour type (all events); and the
    partner preference: the contacts in events of at least MIN_CELLS_FOR_CHANCE cells against the contacts expected
    there by chance, the mean over all ways of shuffling the cell types among the cells of each event (the cells of
    every event and its number of contacts are kept). Observed / expected > 1: within the events, these types touch
    more often than random pairing of the cells would give.

    Events of two cells are left out of that comparison: two touching cells can only form the pair they form, so
    there observed = expected by construction, and they would pull every ratio towards 1. Events without contacts
    take no part (nor do the events left out before the export): this is a preference among the selected multiplets,
    not a binding frequency.

    Expected, exactly: in an event of n cells with c_A cells of type A, a touching pair is of types A and B with
    probability 2 c_A c_B / (n (n - 1)) (A != B), counted once from each side; A–A pairs with probability
    c_A (c_A - 1) / (n (n - 1)), counted twice."""
    types = type_label(cells).to_numpy()
    names = sorted(set(types))
    index = {t: k for k, t in enumerate(names)}
    k = len(names)
    type_idx = np.array([index[t] for t in types], dtype=np.int64)
    row_of = {int(lab): r for r, lab in enumerate(cells["label"].astype(int))}
    frames = cells["frame"].to_numpy(dtype=np.int64)
    observed = np.zeros((k, k))
    observed_3 = np.zeros((k, k))  # in events of at least MIN_CELLS_FOR_CHANCE cells
    event_size = np.bincount(frames, minlength=int(frames.max(initial=-1)) + 1)
    touching: list[list[set[int]]] = [[set() for _ in range(k)] for _ in range(k)]
    pairs_per_event: dict[int, int] = {}
    for a, b in pairs:
        ra, rb = row_of[int(a)], row_of[int(b)]
        ta, tb = type_idx[ra], type_idx[rb]
        observed[ta, tb] += 1
        observed[tb, ta] += 1
        if event_size[frames[ra]] >= MIN_CELLS_FOR_CHANCE:
            observed_3[ta, tb] += 1
            observed_3[tb, ta] += 1
        touching[ta][tb].add(ra)
        touching[tb][ta].add(rb)
        pairs_per_event[int(frames[ra])] = pairs_per_event.get(int(frames[ra]), 0) + 1
    expected = np.zeros((k, k))
    for event, n_pairs in pairs_per_event.items():
        c = np.bincount(type_idx[frames == event], minlength=k).astype(float)
        n = c.sum()
        if n < MIN_CELLS_FOR_CHANCE:
            continue
        p = np.outer(c, c) - np.diag(c)  # c_A c_B, and c_A (c_A - 1) on the diagonal
        expected += n_pairs * 2 * p / (n * (n - 1))
    n_cells = np.bincount(type_idx, minlength=k)
    rows = []
    for i in range(k):
        for j in range(k):
            if observed[i, j] or expected[i, j] > 1e-9:
                rows.append({"cell type": names[i], "neighbour type": names[j], "cells": int(n_cells[i]),
                             "contacts": int(observed[i, j]), "cells touching": len(touching[i][j]),
                             "fraction touching": len(touching[i][j]) / n_cells[i],
                             "contacts (3+ cells)": int(observed_3[i, j]), "expected (3+ cells)": expected[i, j],
                             "observed / expected": (observed_3[i, j] / expected[i, j] if expected[i, j] > 0
                                                     else np.nan)})
    columns = ["cell type", "neighbour type", "cells", "contacts", "cells touching", "fraction touching",
               "contacts (3+ cells)", "expected (3+ cells)", "observed / expected"]
    return pd.DataFrame(rows, columns=columns)


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
