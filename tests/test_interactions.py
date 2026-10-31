"""Touching cells: pairs within the interface reach, neighbour lists, contact table, multiplets."""

import numpy as np
import pandas as pd
import pytest

from imagestream_analysis.interactions import add_neighbours, contact_table, multiplets, touching_pairs
from imagestream_analysis.pipeline import CONTACTS_TABLE, MULTIPLETS_TABLE, run
from test_pipeline import settings_for


@pytest.mark.parametrize("gap, reach, touching", [(0, 1, True), (1, 1, False), (1, 2, True), (2, 2, False)])
def test_reach_is_the_interface_reach(gap, reach, touching):
    """A gap of g background pixels is g + 1 steps: as imagej.dilate, which makes the interfaces."""
    labels = np.zeros((1, 10, 20), dtype=np.int32)
    labels[0, 2:8, 2:8] = 1
    labels[0, 2:8, 8 + gap : 14 + gap] = 2
    assert (touching_pairs(labels, reach).tolist() == [[1, 2]]) is touching


def test_cells_in_other_events_never_touch():
    labels = np.zeros((2, 10, 10), dtype=np.int32)
    labels[0, 2:8, 2:8] = 1
    labels[1, 2:8, 2:8] = 2  # same place, next event
    assert touching_pairs(labels, 3).size == 0


def test_neighbours_contacts_and_multiplets():
    cells = pd.DataFrame({"label": [1, 2, 3, 4, 5], "frame": [0, 0, 0, 1, 2],
                          "Cell type": ["Tumor", "Immune", "Immune", "Tumor", "Artefact"],
                          "Subtype": ["", "A+", "negative", "", ""]})
    pairs = np.array([[1, 2], [1, 3], [2, 3]])
    out = add_neighbours(cells, pairs)
    assert out["Neighbours"].tolist() == [2, 2, 2, 0, 0]
    assert out.loc[0, "Touching (labels)"] == "2,3"
    assert out.loc[0, "Neighbour types"] == "Immune A+; Immune negative"
    table = contact_table(cells, pairs).set_index(["cell type", "neighbour type"])
    assert table.loc[("Tumor", "Immune A+"), "contacts"] == 1
    assert table.loc[("Tumor", "Immune A+"), "cells"] == 2  # two tumor cells, one touches an A+ cell
    assert table.loc[("Tumor", "Immune A+"), "fraction touching"] == 0.5
    assert table.loc[("Immune A+", "Immune negative"), "contacts"] == 1
    m = multiplets(cells, pairs, 4)
    assert m["events with n cells"].tolist() == [1, 2, 0, 1]  # event 3 is empty; events 1, 2 one cell; event 0 three
    assert m["events with n interactions"].tolist() == [3, 0, 0, 1]


def test_pipeline_writes_the_interactions(synthetic):
    s = settings_for(synthetic)
    result = run(s)
    assert len(result.pairs) == 10  # one tumor + one immune cell touching, in every event
    assert (result.events["Interactions"] == 1).all() and (result.events["Artefacts"] == 0).all()
    contacts = pd.read_csv(result.out_dir / f"{CONTACTS_TABLE}.tsv", sep="\t")
    row = contacts[(contacts["cell type"] == "Tumor") & (contacts["neighbour type"] == "Immune")].iloc[0]
    assert row["contacts"] == 10 and row["fraction touching"] == 1
    m = pd.read_csv(result.out_dir / f"{MULTIPLETS_TABLE}.tsv", sep="\t")
    assert m.loc[m["n"] == 2, "events with n cells"].item() == 10


def test_expected_contacts_keep_the_total_and_follow_the_composition():
    """Event 0: tumor 1 touches immune 2; event 1: tumor 3, tumor 4, immune 5, with 3–4 and 4–5 touching."""
    cells = pd.DataFrame({"label": [1, 2, 3, 4, 5], "frame": [0, 0, 1, 1, 1],
                          "Cell type": ["Tumor", "Immune", "Tumor", "Tumor", "Immune"], "Subtype": [""] * 5})
    pairs = np.array([[1, 2], [3, 4], [4, 5]])
    t = contact_table(cells, pairs).set_index(["cell type", "neighbour type"])
    assert t["contacts"].sum() == pytest.approx(t["expected"].sum()) == 6  # 3 pairs, each from both sides
    # event 0 (1 pair, one tumor and one immune): certainly tumor–immune; event 1 (2 pairs; 2 tumor, 1 immune):
    # a pair is tumor–immune with probability 2*2*1/(3*2) = 2/3, tumor–tumor with 2*1/(3*2) = 1/3
    assert t.loc[("Tumor", "Immune"), "expected"] == pytest.approx(1 + 2 * 2 / 3)
    assert t.loc[("Tumor", "Tumor"), "expected"] == pytest.approx(2 * 2 * 1 / 3)  # counted from both sides
    assert t.loc[("Tumor", "Immune"), "observed / expected"] == pytest.approx(2 / (1 + 4 / 3))
    assert ("Immune", "Immune") not in t.index  # one immune cell per event: impossible
