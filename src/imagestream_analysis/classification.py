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
from sklearn.mixture import GaussianMixture

from . import imagej
from .config import ChannelSettings, ClassificationSettings

log = logging.getLogger(__name__)
IMMUNE, TUMOR, ARTEFACT = 1, 2, 3  # values in the class map (macro: changeValues(..., 1) for T cells, 2 for tumor)
CLASS_OF = {"Tumor": 0, "Immune": 1, "Artefact": 2}  # the Class column (macro v1.8: 1 = T cell, 0 = tumor cell)
N_CLUSTERS = 2
KMEANS_SEED = 48  # the macro's randomization_seed (the clustering itself is not bit-identical, see SPEC)
MAX_IMMUNE_TYPES = 6  # immune subtypes that the dashboard colours apart (palette slots 1, 4-8)
THRESHOLD_COLUMNS = ["marker", "role", "cells", "method", "threshold", "positive", "of", "how"]


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


TWO_POPULATIONS = "Two populations"
NOISE_SDS = 2.5  # the upper noise limit of the negative cells: their mode plus this many robust SDs (log scale)
MIN_WEIGHT = 0.005  # a second population needs at least this fraction of the cells
MIN_SEPARATION = 2.0  # ... and its mean this many pooled SDs from the first (Ashman's D)


def two_populations(x: np.ndarray, noise_sds: float = NOISE_SDS) -> tuple[float, str]:
    """The threshold on the log scale x between the negative cells and a second, positive population, and how it was
    found. See manual/analysis/classification.md (technical note on the thresholds) for the reasoning.

    1. The negative cells: the highest peak of a kernel density estimate (Silverman bandwidth) is their mode; their
       spread is estimated from its lower side only (1.4826 x the median distance below the mode), which positive
       cells do not reach. Their upper noise limit is mode + noise_sds x that SD.
    2. A mixture of two normal distributions (EM, started at the mode and at the maximum) is fitted to x. If the upper
       component is credible (at least MIN_WEIGHT of the cells, separated by MIN_SEPARATION pooled SDs, its mean above
       the noise limit), the threshold is where the two weighted densities cross: a cell above it is more likely
       positive than negative. It is never below mode + 2 SD.
    3. Otherwise (no second peak, e.g. a positive tail): the noise limit."""
    x = np.sort(x[np.isfinite(x)])
    n = len(x)
    if n < 10:
        return float(imagej.auto_threshold_dark(x, "MaxEntropy")) if n >= 2 else float("inf"), "too few cells"
    iqr = float(np.subtract(*np.percentile(x, [75, 25])))
    spread = min(float(x.std()), iqr / 1.34) if iqr > 0 else float(x.std())
    bw = max(0.9 * spread * n ** -0.2, 1e-6)
    grid = np.linspace(x[0], x[-1], 1024)
    density = np.zeros(len(grid))
    for chunk in range(0, n, 4096):  # a Gaussian kernel per value, summed
        d = (grid[:, None] - x[None, chunk:chunk + 4096]) / bw
        density += np.exp(-0.5 * d * d).sum(axis=1)
    mode = float(grid[int(np.argmax(density))])
    below = mode - x[x <= mode]
    sd = 1.4826 * float(np.median(below)) if len(below) >= 3 else bw
    sd = sd if sd > 0 else bw
    noise = mode + noise_sds * sd
    gm = GaussianMixture(2, random_state=0, means_init=[[mode], [float(x[-1])]]).fit(x[:, None])
    order = np.argsort(gm.means_.ravel())
    mu = gm.means_.ravel()[order]
    sds = np.sqrt(gm.covariances_.ravel()[order])
    w = gm.weights_[order]
    separation = abs(mu[1] - mu[0]) / float(np.sqrt((sds[0] ** 2 + sds[1] ** 2) / 2))
    if w[1] < MIN_WEIGHT or separation < MIN_SEPARATION or mu[1] < noise:
        return noise, f"no second population: noise limit (mode + {noise_sds:g} SD)"
    g = np.linspace(mu[0], mu[1], 2000)
    p0 = w[0] * np.exp(-0.5 * ((g - mu[0]) / sds[0]) ** 2) / sds[0]
    p1 = w[1] * np.exp(-0.5 * ((g - mu[1]) / sds[1]) ** 2) / sds[1]
    if not (p1 >= p0).any():
        return noise, f"no crossing: noise limit (mode + {noise_sds:g} SD)"
    cross = float(g[int(np.argmax(p1 >= p0))])
    if cross < mode + 2 * sd:
        return mode + 2 * sd, f"two populations ({w[1]:.0%} positive), raised to mode + 2 SD"
    return cross, f"two populations ({w[1]:.0%} positive): crossing"


def log_threshold(values: np.ndarray, method: str, noise_sds: float = NOISE_SDS) -> float:
    return log_threshold_how(values, method, noise_sds)[0]


def log_threshold_how(values: np.ndarray, method: str, noise_sds: float = NOISE_SDS) -> tuple[float, str]:
    """A threshold on intensities that are roughly log-normal (2.8, lines 935-991): the auto-threshold of
    log(v - min + 1), transformed back to intensity.

    '<method> strict' (new): the threshold is raised, if needed, to the median plus STRICT_SDS robust SDs
    (1.4826 x MAD) of the values below it, all on the log scale: no closer to the bulk of the negative cells. On
    synthetic log-normal mixtures this roughly halves the false positives of MaxEntropy, for a few more misses."""
    v = values[np.isfinite(values)]
    if len(v) < 2:
        return float("inf"), "too few cells"
    lo = float(v.min())
    x = np.log(v - lo + 1.0)
    if method == TWO_POPULATIONS:
        t, how = two_populations(x, noise_sds)
        return float(np.exp(t) - 1.0 + lo), how
    t = imagej.auto_threshold_dark(x, method.removesuffix(STRICT))
    if method.endswith(STRICT):
        negative = x[x <= t]
        if len(negative) >= 2:
            median = float(np.median(negative))
            t = max(t, median + STRICT_SDS * 1.4826 * float(np.median(np.abs(negative - median))))
    return float(np.exp(t) - 1.0 + lo), method


def kmeans_classes(df: pd.DataFrame, seed: int = KMEANS_SEED) -> np.ndarray:
    """Cluster on log(immune intensity), log(tumor intensity) and area, each scaled to mean 0, sd 1.

    A cell whose mean intensity is 0 or less has no logarithm. This happens when a channel of its event image is
    empty where the cell lies (seen once: a block of zero pixels in one channel of one IDEAS export). Such a cell
    gets half the lowest positive intensity of the cells instead, so it clusters by its other marker, and the log
    names its event."""
    logs = []
    for name in ("Immune intensity", "Tumor intensity"):
        v = df[name].to_numpy(dtype=float)
        bad = ~(v > 0)
        if bad.any():
            positive = v[v > 0]
            floor = float(positive.min()) / 2 if len(positive) else 1.0
            where = (f" (events {', '.join(str(int(f) + 1) for f in df.loc[bad, 'frame'].head(10))})"
                     if "frame" in df else "")
            log.warning("%d cell(s) with %s 0 or less%s: an empty part of that channel? Clustered with %.3g instead",
                        bad.sum(), name.lower(), where, floor)
            v = np.where(bad, floor, v)
        logs.append(np.log(v))
    features = np.column_stack([*logs, df["Area"].to_numpy(dtype=float)])
    features = (features - features.mean(axis=0)) / features.std(axis=0, ddof=1)  # ImageJ: sample sd
    if not np.isfinite(features).all():
        raise ValueError("Non-finite features for the clustering (a cell with zero or negative marker intensity?)")
    km = KMeans(n_clusters=N_CLUSTERS, n_init=10, random_state=seed).fit(features)
    return np.asarray(km.labels_)


def assign_cell_types(df: pd.DataFrame, clusters: np.ndarray) -> pd.Series:
    """The tumor cells are the cluster with the lowest *sum* of 'Marker ratio over Area' (macro lines 405-420;
    the comment there says average, the code sums: issue I4)."""
    ratio = df["Marker ratio over Area"].replace([np.inf, -np.inf], np.nan)  # a tumor intensity of 0: no ratio
    sums = [ratio[clusters == k].sum() for k in range(N_CLUSTERS)]
    tumor_cluster = 0 if sums[0] < sums[1] else 1
    return pd.Series(np.where(clusters == tumor_cluster, "Tumor", "Immune"), index=df.index)


def classify_pool(datasets: list[tuple[pd.DataFrame, ChannelSettings]],
                  settings: ClassificationSettings) -> tuple[list[pd.DataFrame], pd.DataFrame]:
    """Classify the measured cells of one or more datasets together.

    datasets: (measured cells, resolved channel settings) per dataset; the channels may be numbered differently,
    the markers are matched by role and position. Returns the cells of each dataset with the columns 'Cell type'
    (Tumor, Immune or Artefact), 'Class', '<marker>+' (0/1) per functional marker (with the threshold of the immune
    cells), 'Positive markers' (the number of functional markers a detection is positive in, with the thresholds of
    the artefact rule) and 'Subtype', and the table of thresholds.

    The artefact rule (2026-11-01, after the second user test): a first k-means of all detections only finds the
    immune cells; each subtype marker is thresholded on them, and a detection positive in at least
    artefact_min_positive markers by those thresholds is an artefact. The final k-means runs without the artefacts,
    and the subtype thresholds are computed again on its immune cells. (Thresholds over all detections lay at the
    noise of the many negative cells, so a cell that picked up a little light from a touching neighbour became an
    artefact.)

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
    rows = []

    def threshold(m: Marker, population: pd.Series, purpose: str) -> float:
        manual = settings.thresholds.get(m.name)
        values = pool.loc[population, f"_m {m.key}"].to_numpy(dtype=np.float64)
        t, how = ((float(manual), "set by hand") if manual is not None
                  else log_threshold_how(values, settings.threshold_method))
        positive = pool[f"_m {m.key}"] > t
        rows.append({"marker": m.name, "role": "subtype" if m.role == "functional" else m.role, "cells": purpose,
                     "method": "manual" if manual is not None else settings.threshold_method,
                     "threshold": t, "positive": int((positive & population).sum()), "of": int(population.sum()),
                     "how": how})
        return t

    # 1. the artefacts: positive in too many subtype markers, thresholded on the immune cells of a first clustering
    n_min = settings.artefact_min_positive
    if n_min is None:  # all functional markers, when there are at least two
        n_min = len(markers) if len(markers) >= 2 else 0
    pool["Positive markers"] = 0
    if markers and n_min > 0:
        first = assign_cell_types(pool, kmeans_classes(pool)) == "Immune"
        positive_all = [pool[f"_m {m.key}"] > threshold(m, first, "artefact rule") for m in markers]
        pool["Positive markers"] = np.sum(positive_all, axis=0).astype(int)
    artefact = (pool["Positive markers"] >= n_min if n_min > 0 else pd.Series(False, index=pool.index))
    if artefact.any():
        log.info("%d of %d detections are positive in %d or more subtype markers: artefacts", artefact.sum(),
                 len(pool), n_min)
    elif len(markers) == 1 and settings.artefact_min_positive is None:
        log.info("Artefact detection off: it needs two or more subtype markers (with one, a detection bright in "
                 "every channel counts as a positive immune cell)")
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
        pool[f"{m.name}+"] = (pool[f"_m {m.key}"] > threshold(m, immune, "immune cells")).astype(int)
    pool["Subtype"] = ""
    if functional:
        positive = pool[[f"{m.name}+" for m in functional]].to_numpy(dtype=bool)
        subtype = np.array([" ".join(f"{m.name}+" for m, on in zip(functional, row, strict=True) if on)
                            or "negative" for row in positive], dtype=object)
        pool.loc[immune, "Subtype"] = subtype[immune.to_numpy()]
        n_immune = pool.loc[immune, "Subtype"].nunique()
        if n_immune > MAX_IMMUNE_TYPES:
            log.warning("%d immune subtypes: the dashboard colours %d apart, the others in grey", n_immune,
                        MAX_IMMUNE_TYPES)
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
    extra += [c for c in ("Neighbours", "Touching (labels)", "Neighbour types") if c in cells]
    return pd.concat([out, cells[extra].assign(Edge=cells["Edge"].astype(int))], axis=1)
