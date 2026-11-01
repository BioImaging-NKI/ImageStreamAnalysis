"""Which cells touch which, within each event (ClassifyCellTypes 2.8, lines 479-484 and 609-632; replaces the
Summarize_Interactions notebook).

2.8 uses MorphoLibJ's Region Adjacency Graph on the eroded labels. Here two cells interact when one comes within
`reach` pixels (4-connected steps) of the other: the same reach as the interfaces of the membrane analysis, so that
a cell with an interface always has the partner in its list. Cells are named by their label, which is unique in
the dataset (the 'Label' column of the cell table)."""

from __future__ import annotations

import warnings

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


BOOTSTRAP = 1000  # resamples of the events for the 95 % interval of the subset preference
BOOTSTRAP_SEED = 48
MIN_EXPECTED = 5.0  # fewer expected contacts: no interval (a resample of a handful of contacts says nothing)


def subset_preference(cells: pd.DataFrame, pairs: np.ndarray, bootstrap: int = BOOTSTRAP,
                      seed: int = BOOTSTRAP_SEED, within: str = "event") -> pd.DataFrame:
    """The partner preference of the immune subtypes, among the immune cells: the contacts of each subtype with tumor
    cells and with each subtype, against those expected if the subtypes were shuffled among the immune cells, with
    positions, cell types and contacts kept. The question: when an event offers a choice of immune cells, is the one
    that touches the tumor cell (or another immune cell) more often of this subtype than chance? Only immune cells are
    shuffled, and they are of similar size, so the size bias of contact_table does not arise.

    within = "event" (the default): the subtypes are shuffled among the immune cells of each event, and only events
    with at least two immune cells take part (with one, observed = expected by construction). An immune–tumor contact
    in an event with c_S cells of subtype S among its n immune cells is with S with probability c_S / n; an
    immune–immune contact is S–T with probability 2 c_S c_T / (n (n - 1)), S–S c_S (c_S - 1) / (n (n - 1)), counted
    from both sides as in contact_table. within = "dataset": the same with the counts of the whole dataset and all
    events; when nearly every immune cell touches a tumor cell (gated conjugates) this can hardly depart from 1.

    The 95 % interval: the 2.5 and 97.5 percentiles of observed / expected over `bootstrap` resamples of the events
    (with replacement; seeded, so the tables are reproducible); none below MIN_EXPECTED expected contacts.
    Artefacts take no part."""
    columns = ["cell type", "neighbour type", "contacts", "expected (subtypes shuffled)",
               "observed / expected (subtypes shuffled)", "95% low", "95% high"]
    labels = type_label(cells)
    immune = (cells["Cell type"] == "Immune").to_numpy()
    subtypes = sorted(set(labels[immune]))
    if len(subtypes) < 2:
        return pd.DataFrame(columns=columns)
    k = len(subtypes)
    sub_idx = np.full(len(cells), -1, dtype=np.int64)
    sub_idx[immune] = [subtypes.index(t) for t in labels[immune]]
    is_tumor = (cells["Cell type"] == "Tumor").to_numpy()
    frames = cells["frame"].to_numpy(dtype=np.int64)
    events, event_of = np.unique(frames, return_inverse=True)
    n_ev = len(events)
    counts = np.zeros((n_ev, k))
    np.add.at(counts, (event_of[immune], sub_idx[immune]), 1)
    with_tumor = np.zeros((n_ev, k))  # contacts of subtype S with a tumor cell, per event
    pair_counts = np.zeros((n_ev, k, k))  # immune–immune contacts, both sides
    row_of = {int(lab): r for r, lab in enumerate(cells["label"].astype(int))}
    for a, b in pairs:
        ra, rb = row_of[int(a)], row_of[int(b)]
        sa, sb, e = sub_idx[ra], sub_idx[rb], event_of[ra]
        if sa >= 0 and sb >= 0:
            pair_counts[e, sa, sb] += 1
            pair_counts[e, sb, sa] += 1
        elif sa >= 0 and is_tumor[rb]:
            with_tumor[e, sa] += 1
        elif sb >= 0 and is_tumor[ra]:
            with_tumor[e, sb] += 1

    # per event: the immune cells, and what the shuffle within the event expects (events of one immune cell: equal)
    n_e = counts.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        exp_t_e = np.where(n_e[:, None] > 0, with_tumor.sum(axis=1)[:, None] * counts / n_e[:, None], 0)
        cc_e = counts[:, :, None] * counts[:, None, :] - counts[:, :, None] * np.eye(k)[None]
        exp_ii_e = np.where((n_e > 1)[:, None, None], (pair_counts.sum(axis=(1, 2)) / 2)[:, None, None] * 2 * cc_e
                            / (n_e * (n_e - 1))[:, None, None], 0)
    choice = n_e >= 2  # within="event": only events with a choice of immune cells

    def ratios(w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(observed, expected), (..., k, k + 1): column 0 the tumor cells, then the subtypes; w: event weights."""
        if within == "event":
            w = w * choice
        obs_t = w @ with_tumor
        obs_ii = np.tensordot(w, pair_counts, axes=(-1, 0))
        if within == "event":
            exp_t = w @ exp_t_e
            exp_ii = np.tensordot(w, exp_ii_e, axes=(-1, 0))
        else:
            c = w @ counts
            n = c.sum(axis=-1, keepdims=True)
            n_it = obs_t.sum(axis=-1, keepdims=True)
            n_ii = obs_ii.sum(axis=(-2, -1))[..., None, None] / 2
            with np.errstate(divide="ignore", invalid="ignore"):
                exp_t = n_it * c / n
                cc = c[..., :, None] * c[..., None, :] - np.einsum("...i,ij->...ij", c, np.eye(k))
                exp_ii = n_ii * 2 * cc / (n * (n - 1))[..., None]
        obs = np.concatenate([obs_t[..., :, None], obs_ii], axis=-1)
        exp = np.concatenate([exp_t[..., :, None], exp_ii], axis=-1)
        return obs, exp

    obs, exp = ratios(np.ones(n_ev))
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(n_ev, np.full(n_ev, 1 / n_ev), size=bootstrap).astype(float)
    b_obs, b_exp = ratios(weights)
    with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # a pair without any expected contact in a resample
        b_ratio = np.where(b_exp > 0, b_obs / b_exp, np.nan)
        low, high = np.nanpercentile(b_ratio, [2.5, 97.5], axis=0)
    names = ["Tumor", *subtypes]
    rows = []
    for i, s_name in enumerate(subtypes):
        for j, n_name in enumerate(names):
            if exp[i, j] > 1e-9 or obs[i, j]:
                rows.append({"cell type": s_name, "neighbour type": n_name, "contacts": int(obs[i, j]),
                             "expected (subtypes shuffled)": exp[i, j],
                             "observed / expected (subtypes shuffled)": (obs[i, j] / exp[i, j] if exp[i, j] > 0
                                                                         else np.nan),
                             "95% low": low[i, j] if exp[i, j] >= MIN_EXPECTED else np.nan,
                             "95% high": high[i, j] if exp[i, j] >= MIN_EXPECTED else np.nan})
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
