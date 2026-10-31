import numpy as np
import pandas as pd
import pytest

from imagestream_analysis import io
from imagestream_analysis.config import CellOfInterest, Settings
from imagestream_analysis.pipeline import CELL_TABLE, DOUBLETS_TABLE, EVENT_TABLE, run
from imagestream_analysis.segmentation import clean_labels


def settings_for(synthetic):
    s = Settings()
    s.input.image = synthetic["image"]
    s.input.output = synthetic["tmp"] / "out"
    s.channels.marker_of_interest, s.channels.tumor, s.channels.t_cell = 2, 3, 4
    s.channels.cell_of_interest = CellOfInterest.TUMOR_CELL
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
    s.channels.cell_of_interest = CellOfInterest.TUMOR_CELL
    path = tmp_path / "s.toml"
    path.write_text(s.to_toml(), encoding="utf-8")
    loaded = Settings.load(path)
    assert loaded == s


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
    assert (ev["T cells"] == 1).all() and (ev["Tumor cells"] == 1).all()
    assert len(result.doublets) == 10
    assert (ev["Tumor membrane ratio"] > 1.5).all()  # enrichment at the interface
    out = result.out_dir
    for name in (CELL_TABLE, EVENT_TABLE, DOUBLETS_TABLE):
        assert (out / f"{name}.tsv").exists()
    table = pd.read_csv(out / f"{CELL_TABLE}.tsv", sep="\t")
    assert list(table.columns[:7]) == ["index", "T-cell intensity", "Tumor intensity", "Marker ratio", "Area",
                                       "Marker ratio over Area", "Class"]
    assert (out / "settings.toml").exists() and (out / "run.json").exists()
