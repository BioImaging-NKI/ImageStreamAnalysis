import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from imagestream_analysis.config import CellOfInterest, MembraneSettings
from imagestream_analysis.membrane import (
    INTERFACE_OFFSET,
    _runs,
    analyze_membranes,
    expression_check,
    membrane_profiles,
    outline,
    select,
    straighten,
)


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


def make_cells(bright_partner=False):
    """Event 1: a tumor cell (label 1) touched on its right by an immune cell (3); the marker of interest is bright
    along the tumor cell's interface (or, bright_partner: evenly on the tumor cell, and bright on the immune cell);
    event 2: a lone tumor cell (2)."""
    shape = (2, 80, 80)
    labels = np.zeros(shape, np.int32)
    labels[0][disk(shape[1:], 40, 35, 15)] = 1
    labels[0][disk(shape[1:], 40, 58, 8) & (labels[0] == 0)] = 3
    labels[1][disk(shape[1:], 40, 40, 15)] = 2
    rng = np.random.default_rng(3)
    moi = 100 + rng.normal(0, 5, shape)
    near = (labels[0] == 1) & ndi.binary_dilation(labels[0] == 3, iterations=3)
    if bright_partner:
        moi[0][labels[0] == 3] += 200
    else:
        moi[0][near] += 40  # about 8 sd of the noise
    cells = pd.DataFrame({"frame": [0, 1, 0], "label": [1, 2, 3], "Cell type": ["Tumor", "Tumor", "Immune"],
                          "Subtype": ["", "", "negative"]})
    return labels, cells, moi


def profiles(labels, cells, moi, background=50.0, inside_only=False):
    return membrane_profiles(labels, cells, moi, background, 5, 3, inside_only)


def test_neighbour_tracks():
    labels, cells, moi = make_cells()
    p = profiles(labels, cells, moi)
    assert p.labels.tolist() == [1, 2, 3]
    track = dict(zip(p.labels.tolist(), p.tracks, strict=True))
    assert set(np.unique(track[1])) == {0, 3} and set(np.unique(track[3])) == {0, 1}
    assert (track[2] == 0).all()  # alone in its event
    assert 5 < (track[1] == 3).sum() < 30


def test_enrichment_alignment_and_control():
    labels, cells, moi = make_cells()
    settings = MembraneSettings(control_repeats=10)
    # background 50: the membrane is at 50 above it, the interface at 90
    res = analyze_membranes(cells, profiles(labels, cells, moi), CellOfInterest.TUMOR_CELL, settings)
    t = res.table
    assert len(t) == 2 and res.has_interface.tolist() == [True, False]
    assert 4 < t.loc[0, "stdDevs above mean"] < 30 and t.loc[0, "Ratio interface/non-interface"] > 1.4
    assert np.isnan(t.loc[1, "stdDevs above mean"])
    # the lone cell gets control windows as long as the median interface; its control is about 0
    assert res.window == t.loc[0, "interface length"]
    assert len(res.control_scores) == 10 and abs(t.loc[1, "control stdDevs above mean"]) < 1.5
    # the interface is centred at INTERFACE_OFFSET in the aligned profile
    aligned = res.profiles[0, : res.lengths[0]]
    assert abs(int(np.argmax(np.convolve(aligned, np.ones(5) / 5, "same"))) - INTERFACE_OFFSET) <= 3
    assert res.average["profile"].isna().sum() == 0 or len(res.average) == 0
    assert res.histogram_interface["counts"].sum() == 1 and res.histogram_control["counts"].sum() == 10


def test_any_pair_of_cell_types():
    labels, cells, moi = make_cells()
    p = profiles(labels, cells, moi)
    other_way = analyze_membranes(cells, p, CellOfInterest.IMMUNE_CELL, MembraneSettings())  # immune vs tumor
    assert other_way.table["Label"].tolist() == [3] and other_way.has_interface.tolist() == [True]
    by_subtype = MembraneSettings(cells=["Tumor"], partners=["Immune negative"])
    assert analyze_membranes(cells, p, CellOfInterest.TUMOR_CELL, by_subtype).has_interface.tolist() == [True, False]
    nobody = MembraneSettings(cells=["Tumor"], partners=["Artefact"])
    assert not analyze_membranes(cells, p, CellOfInterest.TUMOR_CELL, nobody).has_interface.any()
    with pytest.raises(ValueError, match="Unknown cell selection"):
        select(cells, ["Macrophage"])


def test_control_is_not_seeded():
    labels, cells, moi = make_cells()
    p = profiles(labels, cells, moi, 100.0)
    a = analyze_membranes(cells, p, CellOfInterest.TUMOR_CELL, MembraneSettings()).control_scores
    b = analyze_membranes(cells, p, CellOfInterest.TUMOR_CELL, MembraneSettings()).control_scores
    assert not np.array_equal(a, b)


def test_a_bright_partner_looks_like_enrichment_and_is_caught():
    """The marker on the immune cell, the tumor cell evenly stained: half the straightened band lies in the
    partner, so the tumor cell's interface looks enriched. The check warns; inside_only removes it."""
    labels, cells, moi = make_cells(bright_partner=True)
    cells["Mean ch2 (MOI)"] = [50.0, 50.0, 250.0]
    s = MembraneSettings(cells=["Tumor"], partners=["Immune"])
    p = profiles(labels, cells, moi)
    fooled = analyze_membranes(cells, p, CellOfInterest.TUMOR_CELL, s).table.loc[0, "stdDevs above mean"]
    inside = profiles(labels, cells, moi, inside_only=True)
    honest = analyze_membranes(cells, inside, CellOfInterest.TUMOR_CELL, s).table.loc[0, "stdDevs above mean"]
    assert fooled > 5 and abs(honest) < 2
    check = expression_check(cells, p, "Mean ch2 (MOI)", ["Tumor"], ["Immune"])
    assert check["check"].startswith("WARNING") and check["partners / measured"] > 1
    assert expression_check(cells, None, "Mean ch2 (MOI)", ["Immune"], ["Tumor"])["check"].startswith("ok")
    # on both at a similar level: membrane doubling, not 'brighter on the partners'
    cells["Mean ch2 (MOI)"] = [100.0, 100.0, 95.0]
    assert "similar level" in expression_check(cells, None, "Mean ch2 (MOI)", ["Tumor"], ["Immune"])["check"]
    # measured cells that hardly carry the marker: not informative (unless the partners are brighter)
    cells["Mean ch2 (MOI)"] = [100.0, 100.0, 4.0]
    assert expression_check(cells, None, "Mean ch2 (MOI)", ["Immune"], ["Immune"])["check"].startswith("NOT INF")
    assert expression_check(cells, None, "Mean ch2 (MOI)", ["Immune"], ["Tumor"])["check"].startswith("WARNING")
