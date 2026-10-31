import numpy as np
import pandas as pd
import pytest

from imagestream_analysis.config import CellOfInterest, MembraneSettings
from imagestream_analysis.membrane import INTERFACE_OFFSET, _runs, analyze_membranes, outline, straighten


def disk(shape, cy, cx, r):
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def test_outline_of_a_disk_has_the_perimeter_and_1_px_steps():
    line = outline(disk((60, 60), 30, 30, 15))
    assert len(line) == pytest.approx(2 * np.pi * 15, rel=0.03)  # along the pixel edges, corners cut
    steps = np.linalg.norm(np.diff(line, axis=0), axis=1)
    assert np.allclose(steps, steps.mean(), rtol=0.05) and steps.mean() == pytest.approx(1, rel=0.03)


def test_straighten_samples_across_the_line():
    mask = disk((60, 60), 30, 30, 15)
    image = mask.astype(float)
    band = straighten(image, outline(mask), 5)
    assert band.shape[0] == 5
    # across the edge: one side inside the cell (1), the other side outside (0)
    assert set(np.round(band[[0, -1]].mean(axis=1))) == {0.0, 1.0}


def test_runs_wrap_around():
    f = np.array([1, 1, 0, 0, 1, 0, 1, 1], bool)
    assert sorted(_runs(f)) == [(4, 1), (6, 4)]
    assert _runs(np.zeros(5, bool)) == []
    assert _runs(np.ones(5, bool)) == [(0, 5)]


def make_cells():
    """Event 1: a tumor cell with an interface of bright MOI on its right side; event 2: a lone tumor cell."""
    shape = (2, 80, 80)
    labels = np.zeros(shape, np.int32)
    labels[0][disk(shape[1:], 40, 35, 15)] = 1
    labels[1][disk(shape[1:], 40, 40, 15)] = 2
    interface = np.zeros(shape, bool)
    interface[0, 30:51, 46:51] = True
    interface &= labels > 0
    rng = np.random.default_rng(3)
    moi = 100 + rng.normal(0, 5, shape)
    moi[0][interface[0]] += 40  # about 8 sd of the noise
    cells = pd.DataFrame({"frame": [0, 1], "label": [1, 2], "Cell type": ["Tumor", "Tumor"]})
    return labels, cells, interface, moi


def test_enrichment_alignment_and_control():
    labels, cells, interface, moi = make_cells()
    settings = MembraneSettings(control_repeats=10)
    # background 50: the membrane is at 50 above it, the interface at 90
    res = analyze_membranes(labels, cells, interface, moi, 50.0, CellOfInterest.TUMOR_CELL, settings)
    t = res.table
    assert len(t) == 2 and res.has_interface.tolist() == [True, False]
    assert 4 < t.loc[0, "stdDevs above mean"] < 20 and t.loc[0, "Ratio interface/non-interface"] > 1.4
    assert np.isnan(t.loc[1, "stdDevs above mean"])
    # the lone cell gets control windows as long as the median interface; its control is about 0
    assert res.window == t.loc[0, "interface length"]
    assert len(res.control_scores) == 10 and abs(t.loc[1, "control stdDevs above mean"]) < 1.5
    # the interface is centred at INTERFACE_OFFSET in the aligned profile
    aligned = res.profiles[0, : res.lengths[0]]
    assert abs(int(np.argmax(np.convolve(aligned, np.ones(5) / 5, "same"))) - INTERFACE_OFFSET) <= 3
    assert res.average["profile"].isna().sum() == 0 or len(res.average) == 0
    assert res.histogram_interface["counts"].sum() == 1 and res.histogram_control["counts"].sum() == 10


def test_control_is_not_seeded():
    labels, cells, interface, moi = make_cells()
    s = MembraneSettings()
    a = analyze_membranes(labels, cells, interface, moi, 100.0, CellOfInterest.TUMOR_CELL, s).control_scores
    b = analyze_membranes(labels, cells, interface, moi, 100.0, CellOfInterest.TUMOR_CELL, s).control_scores
    assert not np.array_equal(a, b)
