import json

import numpy as np
import pandas as pd
import pytest

from imagestream_analysis import io
from imagestream_analysis.config import Settings
from imagestream_analysis.pipeline import CELL_TABLE, DOUBLETS_TABLE, EVENT_TABLE, run
from imagestream_analysis.segmentation import clean_labels


def settings_for(synthetic):
    s = Settings()
    s.input.image = synthetic["image"]
    s.input.output = synthetic["tmp"] / "out"
    s.channels.marker_of_interest, s.channels.tumor, s.channels.immune = [2], 3, [4]
    s.segmentation.method = "labelmap"
    s.segmentation.labelmap_file = synthetic["labelmap"]
    s.membrane.thickness = 2
    s.membrane.analysis = True
    return s


def test_read_stack(synthetic):
    stack = io.read_stack(synthetic["image"], max_events=4)
    assert stack.data.shape == (4, 4, 60, 60)
    assert stack.n_events_in_file == 10
    assert stack.channel_names == ["BF", "MOI", "Tumor", "T cell"]
    assert stack.medians == [800, 10, 20, 20]
    assert stack.pixel_size == pytest.approx(0.5)
    assert np.array_equal(stack.data, synthetic["data"][:4])


def test_montage_one_row_as_the_macro():
    stack = np.arange(3 * 4 * 5).reshape(3, 4, 5)
    montage = io.stack_to_montage(stack)
    assert montage.shape == (4, 15)
    assert np.array_equal(montage[:, 5:10], stack[1])
    assert np.array_equal(io.montage_to_stack(montage, 4, 5), stack)


def test_montage_grid_with_channels_and_empty_tiles():
    stack = np.random.default_rng(1).random((7, 4, 5, 2))
    montage = io.stack_to_montage(stack, columns=3)
    assert montage.shape == (12, 15, 2)  # 3 rows of 3; the last 2 tiles empty
    assert np.array_equal(montage[4:8, 0:5], stack[3])  # row by row
    assert (montage[8:12, 5:15] == 0).all()
    back = io.montage_to_stack(montage, 4, 5)
    assert back.shape == (9, 4, 5, 2) and np.array_equal(back[:7], stack)


def test_grid_columns_about_square():
    assert io.grid_columns(400, 150, 100) == 24  # 24 x 17 tiles: 2400 x 2550 px
    assert io.grid_columns(1, 150, 100) == 1


def test_read_labelmap_grid_montage(tmp_path):
    import tifffile

    labels = np.arange(1, 6)[:, None, None] * np.ones((5, 4, 3), int)
    path = tmp_path / "grid.tif"
    tifffile.imwrite(path, io.stack_to_montage(labels, columns=2).astype(np.uint16))
    assert np.array_equal(io.read_labelmap(path, 5, 4, 3), labels)


def test_cell_across_two_events_is_split():
    labels = np.zeros((3, 6, 6), int)
    labels[0, 2:5, 3:6] = 4  # one Cellpose cell over the border of events 1 and 2
    labels[1, 2:5, 0:3] = 4
    labels[1, 0:2, 4:6] = 2
    labels[2, 1:4, 1:4] = 9
    out = clean_labels(labels, 0, 1, 1000)
    assert out.max() == 4
    assert len(np.unique(out[0])) == 2 and len(np.unique(out[1])) == 3
    assert out[0, 3, 4] == 1 and out[1, 0, 4] == 2 and out[1, 3, 1] == 3 and out[2, 2, 2] == 4


def test_clean_labels_size_filter_only_with_erosion():
    labels = np.zeros((1, 20, 20), int)
    labels[0, 2:5, 2:5] = 5  # 9 px, 1 px after erosion
    labels[0, 8:18, 8:18] = 9  # 100 px, 64 after erosion
    assert np.unique(clean_labels(labels, 0, 50, 1000)).tolist() == [0, 1, 2]
    out = clean_labels(labels, 1, 50, 1000)
    assert np.unique(out).tolist() == [0, 1] and (out == 1).sum() == 64


def test_settings_round_trip(tmp_path):
    s = Settings()
    s.input.image = tmp_path / "a b" / "x.tif"
    s.input.output = tmp_path / "out"
    s.input.max_events = 7
    s.channels.marker_of_interest = [2, "Other"]
    path = tmp_path / "s.toml"
    path.write_text(s.to_toml(), encoding="utf-8")
    loaded = Settings.load(path)
    assert loaded == s
    assert loaded.membrane.max_histogram_value is None and loaded.membrane.histogram_bin_width is None  # automatic
    s.membrane.max_histogram_value, s.membrane.histogram_bin_width = 15.0, 0.5  # explicit values still work
    path.write_text(s.to_toml(), encoding="utf-8")
    assert Settings.load(path) == s
    path.write_text("[membrane]\nmax_histogram_value = 20\nhistogram_bin_width = 1\n", encoding="utf-8")  # ints too
    assert Settings.load(path).membrane.max_histogram_value == 20.0


def test_relative_paths_are_relative_to_the_settings_file(tmp_path):
    path = tmp_path / "s.toml"
    path.write_text("[input]\nimage = 'data/x.tif'\noutput = 'results'\n", encoding="utf-8")
    s = Settings.load(path)
    assert s.input.image == tmp_path / "data" / "x.tif"
    assert s.input.output == tmp_path / "results"


def test_unknown_setting_is_an_error(tmp_path):
    path = tmp_path / "s.toml"
    path.write_text("[input]\nmax_cells = 3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="max_cells"):
        Settings.load(path)


def test_pipeline_on_synthetic_doublets(synthetic):
    result = run(settings_for(synthetic))
    cells = result.cells
    assert len(cells) == 20
    # the large, tumor-marker-bright cell of every event is the tumor cell
    tumor_labels = set(range(1, 21, 2))
    assert set(cells.loc[cells["Cell type"] == "Tumor", "label"]) == tumor_labels
    ev = result.events
    assert (ev["Immune cells"] == 1).all() and (ev["Tumor cells"] == 1).all()
    assert len(result.doublets) == 10
    assert (ev["Tumor membrane ratio (MOI)"] > 1.5).all()  # enrichment at the interface
    out = result.out_dir
    for name in (CELL_TABLE, EVENT_TABLE, DOUBLETS_TABLE):
        assert (out / f"{name}.tsv").exists()
    table = pd.read_csv(out / f"{CELL_TABLE}.tsv", sep="\t")
    assert list(table.columns[:7]) == ["index", "Immune intensity", "Tumor intensity", "Marker ratio", "Area",
                                       "Marker ratio over Area", "Class"]
    assert (out / "settings.toml").exists() and (out / "run.json").exists()


def test_old_settings_files_still_load(tmp_path):
    path = tmp_path / "s.toml"
    path.write_text("[channels]\nt_cell = 6\ncell_of_interest = 'T cell'\n", encoding="utf-8")
    ch = Settings.load(path).channels
    assert ch.immune == [6] and not hasattr(ch, "cell_of_interest")  # dropped (2026-11-03), with a log line
    # one marker of interest (before 2026-11-03): a number or a name becomes a list, "" none
    for old, new in (("2", [2]), ("'CD3'", ["CD3"]), ("''", []), ("[2, 'CD3']", [2, "CD3"])):
        path.write_text(f"[channels]\nmarker_of_interest = {old}\ncell_of_interest = 'Tumor cell'\n",
                        encoding="utf-8")
        assert Settings.load(path).channels.marker_of_interest == new
    path.write_text("[channels]\nt_cell = 6\nimmune = [6]\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not both"):
        Settings.load(path)


def test_channels_by_name_give_the_same_result(synthetic):
    by_number = run(settings_for(synthetic)).cells
    s = settings_for(synthetic)
    s.channels.brightfield, s.channels.marker_of_interest, s.channels.tumor, s.channels.immune = (
        "BF", ["MOI"], "Tumor", ["T cell"])
    by_name = run(s)
    assert by_name.cells.equals(by_number)
    saved = Settings.load(by_name.out_dir / "settings.toml").channels  # saved with the numbers
    assert (saved.brightfield, saved.tumor, saved.immune) == (1, 3, [4])
    s.channels.tumor = "Tumour"
    with pytest.raises(ValueError, match="'Tumour': no channel with this name"):
        run(s)


def test_several_immune_markers_are_summed(synthetic):
    """With the tumor marker also counted as an immune marker, every cell has immune signal: the tumor cells now
    have the larger immune intensity (their own tumor marker), as the sum of the normalized markers."""
    s = settings_for(synthetic)
    one = run(s).cells
    s.channels.immune = [4, 3]
    two = run(s).cells
    assert np.allclose(two["Immune intensity"], one["Immune intensity"] + one["Tumor intensity"], rtol=1e-5)


def test_cellpose_result_is_reused_when_nothing_changed(synthetic, monkeypatch):
    """A second run with the same input and Cellpose settings reuses labelmap_cells.tif; another diameter, or
    reuse_labels = false, runs Cellpose again."""
    from imagestream_analysis import segmentation

    calls = []

    def fake_cellpose(image, channels, settings):
        calls.append(settings.segmentation.diameter)
        return io.stack_to_montage(synthetic["labels"], len(synthetic["labels"])), {"version": "fake"}

    monkeypatch.setattr(segmentation, "run_cellpose", fake_cellpose)
    s = settings_for(synthetic)
    s.segmentation.method, s.segmentation.labelmap_file = "cellpose", None
    s.segmentation.montage_columns = len(synthetic["labels"])
    first = run(s)
    assert calls == [30.0] and (first.out_dir / segmentation.LABELS_INFO).exists()
    second = run(s)
    assert calls == [30.0] and second.cells.attrs["cellpose"]["reused"]
    assert second.cells["label"].tolist() == first.cells["label"].tolist()
    s.segmentation.diameter = 25
    run(s)
    assert calls == [30.0, 25]
    s.segmentation.reuse_labels = False
    run(s)
    assert calls == [30.0, 25, 25]


def test_without_a_marker_of_interest(synthetic):
    """Cell types and contacts only: no marker columns, no marker check; the membrane analysis then needs one."""
    from imagestream_analysis.interfaces import MOI_COLUMNS

    s = settings_for(synthetic)
    s.channels.marker_of_interest = []
    s.membrane.analysis = False
    result = run(s)
    assert len(result.cells) == 20 and len(result.doublets) == 10
    assert not any(c.startswith(tuple(MOI_COLUMNS)) for c in result.events.columns) and result.markers == []
    assert (result.out_dir / "dashboard.html").exists()
    s.membrane.analysis = True
    with pytest.raises(ValueError, match="needs a marker of interest"):
        run(s)


def test_widths_in_micrometres():
    from imagestream_analysis.config import membrane_px, straightening_px

    assert [membrane_px(None, px) for px in (1.0, 0.5, 0.33)] == [2, 3, 5]
    assert [straightening_px(None, px) for px in (1.0, 0.5, 0.33)] == [5, 9, 13]
    assert membrane_px(4, 0.5) == 4 and straightening_px(7, 0.5) == 7  # set in pixels: as it is


def test_pair_table_of_every_membrane_pair(synthetic):
    result = run(settings_for(synthetic))
    t = result.pairs_table.set_index(["membrane of", "interface with"])
    row = t.loc[("Tumor", "Immune")]
    assert row["cells with an interface"] == 10 and row["median stdDevs above mean"] > 2
    assert row["marker check"] and (result.out_dir / "Membrane pairs.tsv").exists()


def test_several_markers_of_interest(synthetic_two_markers):
    """Two markers: the first is on the tumor cells, the second on the immune cells. Each gets the pair of the type
    that carries it (from the data), its own files (named after the marker, every marker alike), its columns in the
    event table, and its rows in Membrane pairs.tsv."""
    s = settings_for(synthetic_two_markers)
    s.channels.marker_of_interest = [2, "Second"]
    result = run(s)
    first, second = result.markers
    assert (first.name, first.carried_by, first.cells, first.partners) == ("MOI", "Tumor", ["Tumor"], ["Immune"])
    assert (second.name, second.carried_by, second.cells, second.partners) == ("Second", "Immune", ["Immune"],
                                                                                ["Tumor"])
    assert second.medians["Immune"] > second.medians["Tumor"] and second.suffix == " - Second"
    assert second.profiles is not None and first.profiles is not None
    assert second.profiles.tracks is first.profiles.tracks  # the outlines and neighbour tracks once
    out = result.out_dir
    for name in ("Membrane intensity analysis", "Histogram_interface_stdDevs_above_mean",
                 "Control windows stdDevs above mean"):
        assert (out / f"{name} - MOI.tsv").exists() and (out / f"{name} - Second.tsv").exists()
        assert not (out / f"{name}.tsv").exists()
    assert (out / "membrane_profiles_normalized - Second.tif").exists()
    own = pd.read_csv(out / "Membrane intensity analysis - Second.tsv", sep="\t")
    assert len(own) == 10 and (own["stdDevs above mean"] > 1).all()  # the 10 immune cells, enriched
    pairs = pd.read_csv(out / "Membrane pairs.tsv", sep="\t")
    assert list(pairs.columns[:4]) == ["marker", "marker on", "membrane of", "interface with"]
    assert pairs.groupby("marker")["marker on"].first().to_dict() == {"MOI": "Tumor", "Second": "Immune"}
    events = pd.read_csv(out / "All cell combinations.tsv", sep="\t")
    assert (events["Tumor membrane ratio (MOI)"] > 1.5).all() and (events["Immune membrane ratio (Second)"] > 1.5).all()
    doublets = pd.read_csv(out / "Single Tumor Single Immune cell.tsv", sep="\t")
    assert {"mean tumor interface (MOI)", "mean immune interface (Second)"} <= set(doublets.columns)
    info = json.loads((out / "run.json").read_text(encoding="utf-8"))
    assert [m["marker on"] for m in info["markers of interest"]] == ["Tumor", "Immune"]
    # a pair in the settings applies to every marker
    s.membrane.cells, s.membrane.partners = ["Tumor"], ["Immune"]
    assert [(m.cells, m.carried_by) for m in run(s).markers] == [(["Tumor"], "Tumor"), (["Tumor"], "Immune")]


def test_cellpose_python_accepts_an_environment_folder(tmp_path):
    """A folder (venv, conda, pixi, or a project with .venv) stands for the Python inside it."""
    import sys

    from imagestream_analysis.segmentation import cellpose_python

    assert cellpose_python(None) == sys.executable and cellpose_python("") == sys.executable
    exe = tmp_path / "env" / ".venv" / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    exe.parent.mkdir(parents=True)
    exe.write_bytes(b"")
    assert cellpose_python(tmp_path / "env") == str(exe)
    assert cellpose_python(tmp_path / "env" / ".venv") == str(exe)
    assert cellpose_python(exe) == str(exe)
    with pytest.raises(ValueError, match="no python executable"):
        cellpose_python(tmp_path)


def test_membrane_pair_is_checked_before_the_run(synthetic):
    """A pair that is no selection stops the run before anything is segmented; a valid subset without cells in this
    dataset only warns: no membrane files for that marker, the rest of the dataset is finished."""
    s = settings_for(synthetic)
    s.membrane.cells = ["Immune T cell"]
    with pytest.raises(ValueError, match="'Immune T cell' is not a subset"):
        run(s)
    s.membrane.cells = ["Immune negative"]  # valid, but no cells here (no subtype markers: no subtypes)
    result = run(s)
    assert result.markers[0].membranes is not None and len(result.markers[0].membranes.table) == 0
    out = result.out_dir
    assert not list(out.glob("Membrane intensity analysis*")) and not list(out.glob("membrane_profiles_normalized*"))
    assert (out / "dashboard.html").exists() and (out / "All cell combinations.tsv").exists()
    assert "no measured cells" in (out / "log.txt").read_text(encoding="utf-8")
