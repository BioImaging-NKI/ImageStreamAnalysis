"""The tumor / immune-cell classification of the measured cells (macro lines 371-507)."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

log = logging.getLogger(__name__)
IMMUNE, TUMOR = 1, 2  # values in the class map (macro: changeValues(..., 1) for T cells, 2 for tumor cells)
N_CLUSTERS = 2
KMEANS_SEED = 48  # the macro's randomization_seed (the clustering itself is not bit-identical, see SPEC)


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


def classify(measured: pd.DataFrame) -> pd.DataFrame:
    """The measured cells (measure.measure_cells) with 'Cell type' and 'Class': 1 = immune cell, 0 = tumor cell (in the
    macro the numbers are the arbitrary k-means cluster numbers)."""
    df = measured.copy()
    if len(df) < N_CLUSTERS:
        raise ValueError(f"Only {len(df)} cells: too few to classify")
    df["Cell type"] = assign_cell_types(df, kmeans_classes(df))
    df["Class"] = (df["Cell type"] == "Immune").astype(int)
    log.info("%d cells: %d immune cells, %d tumor cells", len(df), (df["Class"] == 1).sum(), (df["Class"] == 0).sum())
    return df


def class_map(labels: np.ndarray, cells: pd.DataFrame) -> np.ndarray:
    """Labels -> 1 (immune cell) or 2 (tumor cell)."""
    lut = np.zeros(int(labels.max(initial=0)) + 1, dtype=np.uint8)
    lut[cells["label"].to_numpy()] = np.where(cells["Cell type"] == "Immune", IMMUNE, TUMOR)
    out: np.ndarray = lut[labels]
    return out


def classification_output(cells: pd.DataFrame) -> pd.DataFrame:
    """The table as the macro writes it, plus the event, label and cell type, and the measurements that are new
    in the port (area in µm², edge, the mean of every loaded channel)."""
    out = cells[["Immune intensity", "Tumor intensity", "Marker ratio", "Area", "Marker ratio over Area", "Class"]]
    out.insert(0, "index", np.arange(len(cells)))
    out = out.assign(Frame=cells["frame"] + 1, Label=cells["label"], **{"Cell type": cells["Cell type"]})
    extra = ["Area (µm²)", "Edge"] + [c for c in cells.columns if c.startswith("Mean ch")]
    return pd.concat([out, cells[extra].assign(Edge=cells["Edge"].astype(int))], axis=1)
