import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from imagestream_analysis.config import MembraneSettings
from imagestream_analysis.membrane import (
    INTERFACE_OFFSET,
    _nice_step,
    _runs,
    analyze_membranes,
    expression_check,
    histogram_scale,
    membrane_profiles,
    outline,
    select,
    selector_problem,
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
    res = analyze_membranes(cells, profiles(labels, cells, moi), "Tumor", settings)
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
    other_way = analyze_membranes(cells, p, "Immune", MembraneSettings())  # immune vs tumor
    assert other_way.table["Label"].tolist() == [3] and other_way.has_interface.tolist() == [True]
    by_subtype = MembraneSettings(cells=["Tumor"], partners=["Immune negative"])
    assert analyze_membranes(cells, p, "Tumor", by_subtype).has_interface.tolist() == [True, False]
    nobody = MembraneSettings(cells=["Tumor"], partners=["Artefact"])
    assert not analyze_membranes(cells, p, "Tumor", nobody).has_interface.any()
    with pytest.raises(ValueError, match="Unknown cell selection"):
        select(cells, ["Macrophage"])


def test_control_is_not_seeded():
    labels, cells, moi = make_cells()
    p = profiles(labels, cells, moi, 100.0)
    a = analyze_membranes(cells, p, "Tumor", MembraneSettings()).control_scores
    b = analyze_membranes(cells, p, "Tumor", MembraneSettings()).control_scores
    assert not np.array_equal(a, b)


def test_a_bright_partner_looks_like_enrichment_and_is_caught():
    """The marker on the immune cell, the tumor cell evenly stained: half the straightened band lies in the
    partner, so the tumor cell's interface looks enriched. The check warns; inside_only removes it."""
    labels, cells, moi = make_cells(bright_partner=True)
    cells["Mean ch2 (MOI)"] = [50.0, 50.0, 250.0]
    s = MembraneSettings(cells=["Tumor"], partners=["Immune"])
    p = profiles(labels, cells, moi)
    fooled = analyze_membranes(cells, p, "Tumor", s).table.loc[0, "stdDevs above mean"]
    inside = profiles(labels, cells, moi, inside_only=True)
    honest = analyze_membranes(cells, inside, "Tumor", s).table.loc[0, "stdDevs above mean"]
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


def test_nice_step():
    assert [_nice_step(x) for x in (0.05, 0.1, 0.15, 0.22, 0.3, 0.6, 1.0, 1.5, 3, 7, 12, 30)] == \
        [0.1, 0.1, 0.2, 0.25, 0.5, 1, 1, 2.5, 5, 10, 25, 50]


def test_histogram_scale_automatic():
    rng = np.random.default_rng(1)
    control = rng.normal(0, 1, 1000)
    interface = np.concatenate([rng.normal(3, 2, 200), [500.0], [np.nan]])  # one extreme cell, one without a value
    top, width, n = histogram_scale(interface, control)
    p99 = np.percentile(interface[np.isfinite(interface)], 99)
    assert p99 <= top < p99 + width + 1e-9 and top < 50  # the extreme cell does not set the scale
    assert width in (0.1, 0.2, 0.25, 0.5, 1.0) and 24 <= n <= 60
    assert n * width == pytest.approx(top + 5)  # whole bins from -5
    # at least 5 SD; without any value the macro's 20
    assert histogram_scale(np.zeros(5), np.zeros(5))[0] == 5
    assert histogram_scale(np.zeros(0), np.array([np.nan]))[0] == 20


def test_histogram_scale_explicit_values_are_honoured():
    values = np.linspace(-1, 40, 50)
    assert histogram_scale(values, values, 20.0, 0.25) == (20.0, 0.25, 100)  # the macro's defaults: 100 bins
    top, width, n = histogram_scale(values, values, 10.0, None)
    assert top == 10 and width == 0.25 and n == 60
    top, width, n = histogram_scale(values, values, None, 1.0)
    assert width == 1 and top == pytest.approx(np.ceil(np.percentile(values, 99))) and n == top + 5


def test_histogram_tables_follow_the_settings():
    labels, cells, moi = make_cells()
    res = analyze_membranes(cells, profiles(labels, cells, moi), "Tumor",
                            MembraneSettings(control_repeats=10, max_histogram_value=20.0, histogram_bin_width=0.25))
    assert len(res.histogram_interface) == 100 and res.histogram_interface["bin_start"].iloc[0] == -5
    auto = analyze_membranes(cells, profiles(labels, cells, moi), "Tumor",
                             MembraneSettings(control_repeats=10))
    assert len(auto.histogram_interface) <= 60 and auto.histogram_interface["counts"].sum() == 1


def test_selector_problems():
    """The selections of membrane.cells / partners: subsets as the tables name them, with a suggestion otherwise."""
    markers = ["CD1", "Subtype B"]
    for ok in ("Tumor", "Immune", "Immune negative", "Immune CD1+", "Immune CD1+ Subtype B+", "Subtype B+", "CD1+"):
        assert selector_problem(ok, markers) == "", ok
    assert "did you mean 'CD1+' (every CD1-positive immune cell) or 'Immune CD1+'" in selector_problem("Immune CD1",
                                                                                                      markers)
    assert "did you mean 'Immune CD1+ Subtype B+'" in selector_problem("Immune Subtype B+ CD1+", markers)
    assert "after 'Immune '" in selector_problem("Immune T cells", markers)
    assert "is not a selection" in selector_problem("Macrophage", markers)


def test_pair_per_marker_resolution():
    """Each side on its own: the marker's entry in membrane.pairs, else membrane.cells / partners, else automatic."""
    from imagestream_analysis.membrane import pair_of

    s = MembraneSettings(pairs={"M1": {"cells": ["Tumor"], "partners": ["A+"]}, "M2": {"cells": ["A+"]}})
    assert pair_of(s, "Immune", "M1") == (["Tumor"], ["A+"])
    assert pair_of(s, "Immune", "M2") == (["A+"], ["Tumor"])  # partners: automatic, the other type of the carrier
    assert pair_of(s, "Immune", "M3") == (["Immune"], ["Tumor"])  # no entry: automatic
    s.partners = ["Immune negative"]  # the global setting fills a side that an entry leaves out
    assert pair_of(s, "Immune", "M2") == (["A+"], ["Immune negative"])
    assert pair_of(s, "Tumor", "M3") == (["Tumor"], ["Immune negative"])


def test_excess_over_both_membranes_by_simulation():
    """tools/stats/simulate_excess.py: a marker on both cells without enrichment raises the score but not the excess;
    a marker on the measured cell only gives both about the same; an enrichment on the measured side of a marker on
    both cells is found by the excess."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).parents[1] / "tools" / "stats" / "simulate_excess.py"
    spec = importlib.util.spec_from_file_location("simulate_excess", path)
    sim = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sim)
    for projection in ("max", "sum"):  # the maximum across the band (the macro's) and the sum (dashboard switch)
        score_a, excess_a = sim.simulate("a", n=12, projection=projection)
        score_b, excess_b = sim.simulate("b", n=12, projection=projection)
        score_c, excess_c = sim.simulate("c", n=12, projection=projection)
        assert score_a > 4 and abs(excess_a) < 2  # the partner's membrane alone looks like enrichment; the excess not
        assert excess_b == pytest.approx(score_b, rel=0.15)
        assert excess_c > 0.7 * excess_b and excess_c < score_c - 4


def test_sum_across_the_band():
    """Profiles.sums: the sum over the straightened band, next to the maximum; summed() swaps them."""
    labels = np.zeros((1, 40, 40), np.int32)
    yy, xx = np.mgrid[:40, :40]
    labels[0][(yy - 20) ** 2 + (xx - 20) ** 2 <= 100] = 1
    image = np.ones((1, 40, 40)) * 2.0
    cells = pd.DataFrame({"frame": [0], "label": [1], "Cell type": ["Tumor"]})
    prof = membrane_profiles(labels, cells, image, 0.0, 5, 2)
    assert len(prof.sums) == len(prof.profiles) == 1
    assert np.allclose(prof.profiles[0], 2.0) and np.allclose(prof.sums[0], 10.0)  # uniform: max 2, sum 5 x 2
    assert prof.summed().profiles is prof.sums and prof.summed().tracks is prof.tracks


def _halo_events(n=24, inset=2):
    """Pairs of touching cells whose masks reach `inset` px beyond the stained membrane where they are free (as a
    Cellpose outline does), while at the contact the membrane lies on the contact line. The membrane is stained
    evenly everywhere: no enrichment at the interface."""
    h, w = 60, 90
    labels = np.zeros((n, h, w), np.int32)
    moi = np.zeros((n, h, w), np.float64)
    yy, xx = np.mgrid[:h, :w]
    for t in range(n):
        a, b = disk((h, w), 30, 30, 15), disk((h, w), 30, 56, 15)
        labels[t][a & (xx < 43)] = 1
        labels[t][b & (xx >= 43)] = 2
        # the stained membrane, a ridge of one shape everywhere: `inset` inside the free outline, but on the contact
        # line where the cells touch
        r = 15 - inset
        d = np.minimum(np.abs(np.hypot(yy - 30, xx - 30) - r), np.abs(np.hypot(yy - 30, xx - 56) - r))
        d = np.where((xx >= 39) & (xx <= 47), np.inf, d)
        d = np.minimum(d, np.where(labels[t] > 0, np.abs(xx - 42.5), np.inf))
        moi[t] = 100.0 * np.exp(-d**2 / (2 * 1.2**2))
    cells = pd.DataFrame({"frame": np.repeat(np.arange(n), 2), "label": np.tile([1, 2], n),
                          "Cell type": np.tile(["Tumor", "Immune"], n)})
    return labels, cells, moi


def test_the_band_follows_the_membrane_not_the_outline():
    """A band on the outline, narrower than the gap between a free outline and its membrane, makes an evenly
    stained interface look enriched; centred where the membrane is found (per cell type, free and at contacts), it
    does not. The location reports the membrane about `inset` px inside the free outline and on the contact line."""
    from imagestream_analysis.membrane import location_notes, membrane_profiles_of

    labels, cells, moi = _halo_events()

    def ratio(centre, width=3):
        (prof,) = membrane_profiles_of(labels, cells, [moi], [0.0], width, 2, centre=centre, names=["M"])
        r = []
        for p, tr in zip(prof.profiles, prof.tracks, strict=True):
            r.append(p[tr > 0].mean() / p[tr == 0].mean())
        return float(np.median(r)), prof.location

    on_outline, loc = ratio("outline")
    on_membrane, _ = ratio("membrane")
    assert on_outline > 1.2  # faked enrichment
    assert on_membrane == pytest.approx(1.0, abs=0.1)
    assert ratio("membrane", 9)[0] == pytest.approx(on_membrane, abs=0.03)  # and no longer a matter of the width
    assert loc.offset[("Tumor", "")] == pytest.approx(-2.0, abs=0.5)  # inside the free outline
    assert abs(loc.offset[("Tumor", "Immune")]) < 0.75  # on the contact line
    notes = location_notes(loc, 3, 0.5)
    assert any(level == 30 and "misses part of it" in text for level, text in notes)  # outline: a warning
    table = loc.table(0.5)
    assert set(table["at"]) == {"free membrane", "contact with Immune", "contact with Tumor"}
    assert not table["band placed by it"].any()  # on the outline: nothing moved
    # a cell type that carries the marker at less than half the level of the brightest type: its band stays put
    (prof,) = membrane_profiles_of(labels, cells, [moi], [0.0], 3, 2, centre="membrane", names=["M"],
                                   levels=[{"Tumor": 10.0, "Immune": 30.0}])
    assert ("Tumor", "") not in prof.location.offset and ("Immune", "") in prof.location.offset
    assert prof.location.source[("Tumor", "")].startswith("no marker of interest is carried")
