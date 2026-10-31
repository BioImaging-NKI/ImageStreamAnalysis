"""The measure stage: all-channel means, area in µm², cells cut by the edge of their event image."""

import numpy as np
import pandas as pd
import pytest
import tifffile

from conftest import make_events, write_ideas_zip
from imagestream_analysis import ideas, io
from imagestream_analysis.measure import edge_distance, edge_labels
from imagestream_analysis.pipeline import CELL_TABLE, run
from test_pipeline import settings_for


def test_extents_from_the_ideas_reader(tmp_path):
    data, _ = make_events(4)
    path = write_ideas_zip(tmp_path / "x.zip", data, [1, 2, 3, 4], [1, 2, 3, 4], shrink=(2,))
    stack = ideas.read_dataset(path, "Sample A", [1, 2, 3, 4])
    assert stack.image_extents().tolist() == [[0, 0, 60, 60], [0, 0, 60, 60], [5, 5, 55, 55], [0, 0, 60, 60]]


def test_extents_found_in_a_padded_stack():
    rng = np.random.default_rng(1)
    data = np.full((2, 2, 20, 30), 7, dtype=np.uint16)  # padding = the background 7 in both channels
    data[0, :, 2:18, 5:25] = rng.integers(0, 100, (2, 16, 20))
    data[1] = rng.integers(0, 100, (2, 20, 30))
    data[1, :, 0, 0] = 7  # a single pixel at the background value is no padding
    stack = io.EventStack(data, ["a", "b"], [7.0, 7.0], 1.0, "x", 2)
    assert stack.image_extents().tolist() == [[2, 5, 18, 25], [0, 0, 20, 30]]
    stack.medians = None  # unknown background: the whole frame
    assert stack.image_extents().tolist() == [[0, 0, 20, 30], [0, 0, 20, 30]]


def test_edge_cells():
    labels = np.zeros((1, 20, 20), dtype=np.int32)
    labels[0, 1:4, 8:12] = 1  # 1 px from the top of the frame
    labels[0, 8:12, 8:12] = 2  # in the middle
    labels[0, 8:12, 15:17] = 3  # touches the right edge of an image that ends at x = 17
    full, cropped = np.array([[0, 0, 20, 20]]), np.array([[0, 0, 20, 17]])
    assert edge_distance(labels, full).tolist() == [1, 8, 3]
    assert edge_labels(labels, full, 1).tolist() == [True, False, False]
    assert edge_labels(labels, full, 0).tolist() == [False, False, False]
    assert edge_labels(labels, cropped, 1).tolist() == [True, False, True]


def test_cell_table_has_all_channels_area_and_edge(synthetic):
    s = settings_for(synthetic)
    s.membrane.analysis = False
    result = run(s)
    table = pd.read_csv(result.out_dir / f"{CELL_TABLE}.tsv", sep="\t")
    names = ["BF", "MOI", "Tumor", "T cell"]
    assert [c for c in table.columns if c.startswith("Mean ch")] == [f"Mean ch{i + 1} ({n})" for i, n in
                                                                      enumerate(names)]
    assert table["Area (µm²)"].to_numpy() == pytest.approx(table["Area"].to_numpy() * 0.25)
    assert (table["Edge"] == 0).all()
    # background-subtracted: a tumor cell has about 400 + 20 t above the background in the tumor channel
    tumor = table[table["Cell type"] == "Tumor"]
    assert (tumor["Mean ch3 (Tumor)"] > 380).all()
    assert tumor["Mean ch4 (T cell)"].abs().max() < 5


def test_edge_cells_are_flagged_or_left_out(synthetic):
    s = settings_for(synthetic)
    s.membrane.analysis = False
    labels = synthetic["labels"].copy()
    labels[0, 0:12, 40:58] = 99  # a cell cut by the top of event 0
    tifffile.imwrite(synthetic["labelmap"], labels.astype(np.uint16), imagej=True, metadata={"axes": "TYX"})
    cells = run(s).cells
    assert len(cells) == 21 and cells["Edge"].sum() == 1
    s.segmentation.exclude_edge_cells = True
    cells = run(s).cells
    assert len(cells) == 20 and not cells["Edge"].any()
