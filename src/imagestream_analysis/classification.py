"""Per-cell measurements and the tumor / T-cell classification (macro lines 325-507)."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from sklearn.cluster import KMeans

log = logging.getLogger(__name__)
T_CELL, TUMOR = 1, 2  # values in the class map (macro: changeValues(..., 1) for T cells, 2 for tumor cells)
N_CLUSTERS = 2
KMEANS_SEED = 48  # the macro's randomization_seed (the clustering itself is not bit-identical, see SPEC)


def cell_table(labels: np.ndarray, tumor: np.ndarray, t_cell: np.ndarray) -> pd.DataFrame:
    """One row per cell, ordered by event and then label (the order of the macro's ROI Manager).

    labels: (T, Y, X) cleaned labels 1..N; tumor, t_cell: normalized marker channels (T, Y, X)."""
    n = int(labels.max(initial=0))
    index = np.arange(1, n + 1)
    area = ndi.sum_labels(np.ones_like(labels, dtype=np.int64), labels, index)
    t_int = ndi.mean(t_cell, labels, index)
    tu_int = ndi.mean(tumor, labels, index)
    frame_of = np.full(n + 1, -1, dtype=np.int64)
    present = labels.reshape(labels.shape[0], -1)
    for t in range(labels.shape[0]):  # first event in which each label occurs
        vals = np.unique(present[t])
        vals = vals[(vals > 0) & (frame_of[vals] < 0)]
        frame_of[vals] = t
    df = pd.DataFrame({"label": index, "frame": frame_of[1:], "T-cell intensity": t_int,
                       "Tumor intensity": tu_int, "Area": area.astype(np.int64)})
    df = df.sort_values(["frame", "label"], kind="stable").reset_index(drop=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["Marker ratio"] = df["T-cell intensity"] / df["Tumor intensity"]
        df["Marker ratio over Area"] = df["Marker ratio"] / df["Area"]
    return df


def kmeans_classes(df: pd.DataFrame, seed: int = KMEANS_SEED) -> np.ndarray:
    """Cluster on log(T-cell intensity), log(tumor intensity) and area, each scaled to mean 0, sd 1."""
    with np.errstate(divide="ignore", invalid="ignore"):
        features = np.column_stack([np.log(df["T-cell intensity"]), np.log(df["Tumor intensity"]), df["Area"]])
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
    return pd.Series(np.where(clusters == tumor_cluster, "Tumor", "T cell"), index=df.index)


def classify(labels: np.ndarray, tumor: np.ndarray, t_cell: np.ndarray) -> pd.DataFrame:
    """The 'Cell classification' table. Class: 1 = T cell, 0 = tumor cell (in the macro the numbers are the
    arbitrary k-means cluster numbers)."""
    df = cell_table(labels, tumor, t_cell)
    if len(df) < N_CLUSTERS:
        raise ValueError(f"Only {len(df)} cells: too few to classify")
    df["Cell type"] = assign_cell_types(df, kmeans_classes(df))
    df["Class"] = (df["Cell type"] == "T cell").astype(int)
    log.info("%d cells: %d T cells, %d tumor cells", len(df), (df["Class"] == 1).sum(), (df["Class"] == 0).sum())
    return df


def class_map(labels: np.ndarray, cells: pd.DataFrame) -> np.ndarray:
    """Labels -> 1 (T cell) or 2 (tumor cell)."""
    lut = np.zeros(int(labels.max(initial=0)) + 1, dtype=np.uint8)
    lut[cells["label"].to_numpy()] = np.where(cells["Cell type"] == "T cell", T_CELL, TUMOR)
    out: np.ndarray = lut[labels]
    return out


def classification_output(cells: pd.DataFrame) -> pd.DataFrame:
    """The table as the macro writes it, plus the event, label and cell type."""
    out = cells[["T-cell intensity", "Tumor intensity", "Marker ratio", "Area", "Marker ratio over Area", "Class"]]
    out.insert(0, "index", np.arange(len(cells)))
    out = out.assign(Frame=cells["frame"] + 1, Label=cells["label"], **{"Cell type": cells["Cell type"]})
    return out
