"""The classification of the measured cells, in two steps (SPEC, "Merging ClassifyCellTypes 2.8"):

1. lineage: tumor or immune cell, by k-means as the macro v1.8 (lines 371-507), after leaving out the artefacts
   (detections positive in too many functional markers, thresholded on all detections);
2. immune subtypes: positive or negative for each functional marker, by a threshold computed on the immune cells
   (ClassifyCellTypes 2.8, lines 929-1000).

Several datasets can be classified together (a pool): one set of thresholds and one clustering for all their cells.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from . import imagej
from .config import ChannelSettings, ClassificationSettings

log = logging.getLogger(__name__)
IMMUNE, TUMOR, ARTEFACT = 1, 2, 3  # values in the class map (macro: changeValues(..., 1) for T cells, 2 for tumor)
CLASS_OF = {"Tumor": 0, "Immune": 1, "Artefact": 2}  # the Class column (macro v1.8: 1 = T cell, 0 = tumor cell)
N_CLUSTERS = 2
KMEANS_SEED = 48  # the macro's randomization_seed (the clustering itself is not bit-identical, see SPEC)
MAX_TYPES = 8  # cell types that the dashboard can colour apart
THRESHOLD_COLUMNS = ["marker", "role", "cells", "method", "threshold", "positive", "of"]


@dataclass
class Marker:
    key: str  # "tumor", "immune 1", "functional 2": the same marker in every dataset of a pool
    role: str  # tumor, immune or functional
    name: str  # channel name (of the first dataset)


def marker_columns(cells: pd.DataFrame, ch: ChannelSettings) -> list[tuple[Marker, str]]:
    """The markers of the classification and their 'Mean ch<N> (<name>)' column in the measured table. ch must be
    resolved (channel numbers)."""
    out = []
    refs = [("tumor", ch.tumor)] + [(f"immune {i + 1}", c) for i, c in enumerate(ch.immune)]
    refs += [(f"functional {i + 1}", c) for i, c in enumerate(ch.functional)]
    for key, number in refs:
        prefix = f"Mean ch{int(number)} ("
        column = next((c for c in cells.columns if c.startswith(prefix)), None)
        if column is None:
            raise ValueError(f"No measurement of channel {number} ({key} marker)")
        out.append((Marker(key, key.split()[0], column[len(prefix) : -1]), column))
    return out


STRICT = " strict"
STRICT_SDS = 2.0  # "<method> strict": at least this many robust SDs above the median of the negative cells


def log_threshold(values: np.ndarray, method: str) -> float:
    """A threshold on intensities that are roughly log-normal (2.8, lines 935-991): the auto-threshold of
    log(v - min + 1), transformed back to intensity.

    '<method> strict' (new): the threshold is raised, if needed, to the median plus STRICT_SDS robust SDs
    (1.4826 x MAD) of the values below it, all on the log scale: no closer to the bulk of the negative cells. On
    synthetic log-normal mixtures this roughly halves the false positives of MaxEntropy, for a few more misses."""
    v = values[np.isfinite(values)]
    if len(v) < 2:
        return float("inf")
    lo = float(v.min())
    x = np.log(v - lo + 1.0)
    t = imagej.auto_threshold_dark(x, method.removesuffix(STRICT))
    if method.endswith(STRICT):
        negative = x[x <= t]
        if len(negative) >= 2:
            median = float(np.median(negative))
            t = max(t, median + STRICT_SDS * 1.4826 * float(np.median(np.abs(negative - median))))
    return float(np.exp(t) - 1.0 + lo)


def kmeans_classes(df: pd.DataFrame, seed: int = KMEANS_SEED) -> np.ndarray:
    """Cluster on log(immune intensity), log(tumor intensity) and area, each scaled to mean 0, sd 1."""
    with np.errstate(divide="ignore", invalid="ignore"):
        features = np.column_stack([np.log(df["Immune intensity"]), np.log(df["Tumor intensity"]), df["Area"]])
    features = (features - features.mean(axis=0)) / features.std(axis=0, ddof=1)  # ImageJ: sample sd
    if not np.isfinite(features).all():
        raise ValueError("Non-finite features for the clustering (a cell with zero or negative marker intensity?)")
    km = KMeans(n_clusters=N_CLUSTERS, n_init=10, random_state=seed).fit(features)
    return np.asarray(km.labels_)


def assign_cell_types(df: pd.DataFrame, clusters: np.ndarray) -> pd.Series:
    """The tumor cells are the cluster with the lowest *sum* of 'Marker ratio over Area' (macro lines 405-420;
    the comment there says average, the code sums: issue I4)."""
    sums = [df.loc[clusters == k, "Marker ratio over Area"].sum() for k in range(N_CLUSTERS)]
    tumor_cluster = 0 if sums[0] < sums[1] else 1
    return pd.Series(np.where(clusters == tumor_cluster, "Tumor", "Immune"), index=df.index)


def classify_pool(datasets: list[tuple[pd.DataFrame, ChannelSettings]],
                  settings: ClassificationSettings) -> tuple[list[pd.DataFrame], pd.DataFrame]:
    """Classify the measured cells of one or more datasets together.

    datasets: (measured cells, resolved channel settings) per dataset; the channels may be numbered differently,
    the markers are matched by role and position. Returns the cells of each dataset with the columns 'Cell type'
    (Tumor, Immune or Artefact), 'Class', '<marker>+' (0/1) per functional marker (with the threshold of the immune
    cells), 'Positive markers' (the number of functional markers a detection is positive in, with the thresholds of
    all detections: the artefact rule) and 'Subtype', and the table of thresholds.

    The tumor and immune markers are not thresholded: on the test data an automatic threshold over all detections
    calls far more cells positive than the k-means finds, so they would make most detections artefacts."""
    if not datasets:
        return [], pd.DataFrame(columns=THRESHOLD_COLUMNS)
    markers = [m for m, _ in marker_columns(*datasets[0]) if m.role == "functional"]
    parts = []
    for k, (cells, ch) in enumerate(datasets):
        own = [mc for mc in marker_columns(cells, ch) if mc[0].role == "functional"]
        values = {f"_m {m.key}": cells[column].to_numpy(dtype=np.float64)
                  for m, (_, column) in zip(markers, own, strict=True)}
        parts.append(cells.assign(_set=k, **values))
    pool = pd.concat(parts, ignore_index=True)
    if len(pool) < N_CLUSTERS:
        raise ValueError(f"Only {len(pool)} cells: too few to classify")
    everyone = pd.Series(True, index=pool.index)
    rows = []

    def threshold(m: Marker, population: pd.Series) -> float:
        manual = settings.thresholds.get(m.name)
        values = pool.loc[population, f"_m {m.key}"].to_numpy(dtype=np.float64)
        t = float(manual) if manual is not None else log_threshold(values, settings.threshold_method)
        positive = pool[f"_m {m.key}"] > t
        rows.append({"marker": m.name, "role": m.role,
                     "cells": "all cells" if population is everyone else "immune cells",
                     "method": "manual" if manual is not None else settings.threshold_method,
                     "threshold": t, "positive": int((positive & population).sum()), "of": int(population.sum())})
        return t

    # 1. the functional markers on all detections: the artefacts are positive in too many of them
    positive_all = [pool[f"_m {m.key}"] > threshold(m, everyone) for m in markers]
    pool["Positive markers"] = np.sum(positive_all, axis=0).astype(int) if markers else 0
    n_min = settings.artefact_min_positive
    if n_min is None:  # all functional markers, when there are at least two
        n_min = len(markers) if len(markers) >= 2 else 0
    artefact = (pool["Positive markers"] >= n_min if n_min > 0 else pd.Series(False, index=pool.index))
    if artefact.any():
        log.info("%d of %d detections are positive in %d or more functional markers: artefacts", artefact.sum(),
                 len(pool), n_min)
    # 2. tumor or immune cell: k-means on the other detections
    cells_ok = pool[~artefact]
    if len(cells_ok) < N_CLUSTERS:
        raise ValueError(f"Only {len(cells_ok)} cells that are not artefacts: too few to classify")
    pool["Cell type"] = "Artefact"
    pool.loc[cells_ok.index, "Cell type"] = assign_cell_types(cells_ok, kmeans_classes(cells_ok))
    pool["Class"] = pool["Cell type"].map(CLASS_OF)
    # 3. immune subtypes: the functional markers, thresholded on the immune cells
    immune = pool["Cell type"] == "Immune"
    functional = markers
    for m in functional:
        pool[f"{m.name}+"] = (pool[f"_m {m.key}"] > threshold(m, immune)).astype(int)
    pool["Subtype"] = ""
    if functional:
        positive = pool[[f"{m.name}+" for m in functional]].to_numpy(dtype=bool)
        subtype = np.array([" ".join(f"{m.name}+" for m, on in zip(functional, row, strict=True) if on)
                            or "negative" for row in positive], dtype=object)
        pool.loc[immune, "Subtype"] = subtype[immune.to_numpy()]
        n_types = pool.loc[immune, "Subtype"].nunique() + 1 + int(artefact.any())
        if n_types > MAX_TYPES:
            log.warning("%d cell types (tumor, immune subtypes, artefacts): the dashboard colours only %d apart",
                        n_types, MAX_TYPES)
    counts = pool["Cell type"].value_counts()
    log.info("%d cells: %d immune cells, %d tumor cells%s", len(pool), counts.get("Immune", 0),
             counts.get("Tumor", 0), f", {counts['Artefact']} artefacts" if "Artefact" in counts else "")
    pool = pool.drop(columns=[c for c in pool.columns if c.startswith("_m ")])
    out = [pool[pool["_set"] == k].drop(columns="_set").reset_index(drop=True) for k in range(len(datasets))]
    return out, pd.DataFrame(rows, columns=THRESHOLD_COLUMNS)


def classify(measured: pd.DataFrame, channels: ChannelSettings,
             settings: ClassificationSettings) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One dataset on its own: classify_pool with one dataset."""
    out, thresholds = classify_pool([(measured, channels)], settings)
    return out[0], thresholds


def class_map(labels: np.ndarray, cells: pd.DataFrame) -> np.ndarray:
    """Labels -> 1 (immune cell), 2 (tumor cell) or 3 (artefact)."""
    lut = np.zeros(int(labels.max(initial=0)) + 1, dtype=np.uint8)
    values = cells["Cell type"].map({"Immune": IMMUNE, "Tumor": TUMOR, "Artefact": ARTEFACT})
    lut[cells["label"].to_numpy()] = values.to_numpy(dtype=np.uint8)
    out: np.ndarray = lut[labels]
    return out


def classification_output(cells: pd.DataFrame) -> pd.DataFrame:
    """The table as the macro writes it, plus the event, label and cell type, and the measurements that are new
    in the port (area in µm², edge, the mean of every loaded channel, the positive markers and the subtype)."""
    out = cells[["Immune intensity", "Tumor intensity", "Marker ratio", "Area", "Marker ratio over Area", "Class"]]
    out.insert(0, "index", np.arange(len(cells)))
    out = out.assign(Frame=cells["frame"] + 1, Label=cells["label"], **{"Cell type": cells["Cell type"]})
    extra = ["Area (µm²)", "Edge"] + [c for c in cells.columns if c.startswith("Mean ch")]
    extra += [c for c in cells.columns if c.endswith("+")] + ["Positive markers", "Subtype"]
    return pd.concat([out, cells[extra].assign(Edge=cells["Edge"].astype(int))], axis=1)
