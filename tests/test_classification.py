"""The two-step classification: artefacts, tumor / immune cells by k-means, immune subtypes by thresholds; pools."""

import numpy as np
import pandas as pd
import pytest

from imagestream_analysis import imagej
from imagestream_analysis.classification import classify, classify_pool, log_threshold
from imagestream_analysis.config import ChannelSettings, ClassificationSettings

# Fiji (2026-10-31): setAutoThreshold("<method> dark") + getThreshold on these values as a 32-bit 1-pixel-high image,
# and on log(v - min + 1) of them, as ClassifyCellTypes 2.8 thresholds its markers. (On 12 larger log-normal mixtures
# the port also equalled Fiji to 1e-11.)
VALUES = [2.54, 5.72, 5.31, 1.65, 1.98, 1.63, 3.63, 2.38, 4.04, 0.14, 6.39, 2.31, 3.88, 2.24, 1.85, 3.39, 4.23, 2.13,
          2.22, 3.9, 1.16, 0.45, 3.25, 1.43, 0.08, 1.24, 1.72, 0.78, 0.47, 2.55, 4.42, 2.08, 1.33, 3.23, 3.97, 1.97,
          3.57, 4.8, 2.13, 1.24, 58.6, 56.81, 73.92, 35.14, 42.77, 40.46, 30.45, 54.71, 61.97, 41.74, 80.74, 67.87]
FIJI = {("raw", "MaxEntropy"): 4.508392038, ("raw", "Otsu"): 6.722588058,
        ("log", "MaxEntropy"): 1.450256393, ("log", "Otsu"): 2.002735018}


@pytest.mark.parametrize("scale, method", list(FIJI))
def test_auto_threshold_equals_fiji(scale, method):
    v = np.array(VALUES, dtype=np.float32)
    if scale == "log":
        v = np.log(v - v.min() + 1).astype(np.float32)
    assert imagej.auto_threshold_dark(v, method) == pytest.approx(FIJI[scale, method], abs=1e-6)


def test_max_entropy_splits_two_peaks():
    """MaxEntropy cuts in the upper tail of the large (negative) peak: above most of it, below the positives."""
    rng = np.random.default_rng(0)
    v = np.concatenate([rng.normal(50, 5, 1000), rng.normal(150, 10, 200)])
    hist = np.bincount(np.clip(v, 0, 255).astype(int), minlength=256)
    assert 55 < imagej.max_entropy_bin(hist) < 120
    assert 55 < imagej.auto_threshold_dark(v, "MaxEntropy") < 120


def test_log_threshold_on_log_normal_populations():
    rng = np.random.default_rng(1)
    v = np.concatenate([rng.lognormal(2, 0.3, 800), rng.lognormal(5, 0.3, 200)]) - 3  # background-subtracted
    t = log_threshold(v, "MaxEntropy")
    assert np.exp(2.0) - 3 < t < np.exp(4.1) - 3  # above the median of the negatives, below the positives


def measured(rng, n_tumor=200, n_immune=200, n_artefact=0, numbers=(3, 4, 5, 6)):
    """A measured cell table: tumor cells large and tumor-marker bright, immune cells small and immune-marker
    bright; among the immune cells 40 % positive for marker A, another 30 % for marker B (never both); artefacts
    bright in A and B. numbers: the channel numbers of tumor, immune, A, B."""
    n = n_tumor + n_immune + n_artefact
    kind = np.array(["T"] * n_tumor + ["I"] * n_immune + ["X"] * n_artefact)
    tumor = np.where(kind == "T", rng.lognormal(1.0, 0.2, n), rng.lognormal(-0.5, 0.2, n))
    immune = np.where(kind == "I", rng.lognormal(1.0, 0.2, n), rng.lognormal(-0.5, 0.2, n))
    area = np.where(kind == "T", rng.normal(600, 60, n), rng.normal(200, 20, n))
    u = rng.random(n)
    a_pos = ((kind == "I") & (u < 0.4)) | (kind == "X")
    b_pos = ((kind == "I") & (u >= 0.4) & (u < 0.7)) | (kind == "X")
    a = np.where(a_pos, rng.lognormal(5, 0.3, n), rng.lognormal(1.5, 0.3, n))
    b = np.where(b_pos, rng.lognormal(5.5, 0.3, n), rng.lognormal(1.5, 0.3, n))
    df = pd.DataFrame({"label": np.arange(1, n + 1), "frame": np.arange(n) // 3, "Immune intensity": immune,
                       "Tumor intensity": tumor, "Area": area})
    df["Marker ratio"] = df["Immune intensity"] / df["Tumor intensity"]
    df["Marker ratio over Area"] = df["Marker ratio"] / df["Area"]
    t, i, ca, cb = numbers
    df[f"Mean ch{t} (Tumor)"] = tumor * 100
    df[f"Mean ch{i} (Immune)"] = immune * 100
    df[f"Mean ch{ca} (A)"] = a
    df[f"Mean ch{cb} (B)"] = b
    channels = ChannelSettings(tumor=t, immune=[i], functional=[ca, cb])
    return df, channels, kind, a_pos, b_pos


def test_lineage_and_subtypes():
    df, ch, kind, a_pos, b_pos = measured(np.random.default_rng(2))
    cells, thresholds = classify(df, ch, ClassificationSettings())
    assert (cells["Cell type"].to_numpy() == np.where(kind == "T", "Tumor", "Immune")).all()
    assert (cells["Class"] == (kind == "I").astype(int)).all()
    immune = kind == "I"
    for marker, pos in (("A+", a_pos), ("B+", b_pos)):  # every positive found; MaxEntropy is liberal on negatives
        assert (cells.loc[immune, marker] >= pos[immune]).all()
        assert (cells.loc[immune, marker] == pos[immune]).mean() > 0.85
    sub = cells.loc[immune, "Subtype"]
    assert {"A+", "B+", "negative"} <= set(sub) and (cells.loc[~immune, "Subtype"] == "").all()
    assert (sub == "A+ B+").sum() == (cells.loc[immune, "A+"] & cells.loc[immune, "B+"]).sum()
    assert list(thresholds["cells"]) == ["all cells", "all cells", "immune cells", "immune cells"]
    assert list(thresholds["marker"]) == ["A", "B", "A", "B"]
    assert thresholds.loc[2, "of"] == immune.sum()


def test_artefacts_are_left_out_of_the_clustering():
    df, ch, kind, _, _ = measured(np.random.default_rng(3), n_artefact=25)
    cells, _ = classify(df, ch, ClassificationSettings())
    assert "Artefact" not in set(cells["Cell type"])  # off by default
    cells, _ = classify(df, ch, ClassificationSettings(artefact_min_positive=2))
    assert (cells.loc[kind == "X", "Cell type"] == "Artefact").all()
    assert (cells["Cell type"] == "Artefact").sum() < 25 + 15  # a few immune cells positive in both by chance
    assert (cells.loc[kind == "X", "Class"] == 2).all() and (cells.loc[kind == "X", "Subtype"] == "").all()
    # MaxEntropy cuts in the tail of the negatives, so a few true cells are positive in both by chance
    assert (cells.loc[kind == "T", "Cell type"] == "Tumor").mean() > 0.95


def test_manual_threshold_overrides_the_automatic_one():
    df, ch, kind, _, _ = measured(np.random.default_rng(4))
    cells, thresholds = classify(df, ch, ClassificationSettings(thresholds={"A": 1e9}))
    assert cells["A+"].sum() == 0 and (thresholds.loc[thresholds["marker"] == "A", "method"] == "manual").all()
    assert cells["B+"].sum() > 0


def test_a_pool_shares_thresholds_and_matches_channels_by_role():
    rng = np.random.default_rng(5)
    first, ch1, *_ = measured(rng)
    second, ch2, *_ = measured(rng, n_tumor=50, n_immune=300, numbers=(9, 2, 11, 7))  # numbered differently
    (c1, c2), thresholds = classify_pool([(first, ch1), (second, ch2)], ClassificationSettings())
    assert len(c1) == len(first) and len(c2) == len(second)
    assert thresholds["of"].iloc[0] == len(first) + len(second)
    t_a = thresholds.query("marker == 'A' and cells == 'immune cells'")["threshold"].item()
    assert ((c2["Mean ch11 (A)"] > t_a) == c2["A+"].astype(bool)).all()  # the pool's threshold in both datasets
    alone = classify(second, ch2, ClassificationSettings())[1]
    assert alone.query("marker == 'A' and cells == 'immune cells'")["threshold"].item() != pytest.approx(t_a)
