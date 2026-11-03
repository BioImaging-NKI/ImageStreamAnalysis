"""Comparison with the output of the Fiji macro v1.8 on a confidential test data set.

Nothing about those data is in this repository: the test reads a gitignored file, `data/reference_test.toml` (or the
file in ISA_REFERENCE), and is skipped without it:

    [reference]
    settings = 'reference.toml'            # run settings for the test data (relative to this file)
    fiji_output = 'D:/path/to/the/Fiji/output/folder'  # relative to this file, or absolute

The run settings must use the labelmap of the Fiji run (segmentation.method = "labelmap") and the same number of
events, so that the cells are the same.
"""

import os
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from imagestream_analysis.config import CellOfInterest, Settings
from imagestream_analysis.pipeline import AVERAGE_PROFILE, CELL_TABLE, DOUBLETS_TABLE, EVENT_TABLE, MEMBRANE_TABLE, run

CONFIG = Path(os.environ.get("ISA_REFERENCE", Path(__file__).parents[1] / "data" / "reference_test.toml"))

pytestmark = [pytest.mark.reference, pytest.mark.skipif(not CONFIG.exists(), reason=f"no {CONFIG}")]


def _config() -> tuple[Settings, Path]:
    with CONFIG.open("rb") as f:
        ref = tomllib.load(f)["reference"]
    return Settings.load(CONFIG.parent / ref["settings"]), CONFIG.parent / ref["fiji_output"]


@pytest.fixture(scope="module")
def reference() -> Path:
    return _config()[1]


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    settings, _ = _config()
    settings.input.output = tmp_path_factory.mktemp("reference")
    settings.membrane.analysis = True
    return run(settings)


# the port says immune cell where the macro v1.8 says T cell
RENAMED = {"T-cell intensity": "Immune intensity", "T cells": "Immune cells", "mean Tcell interface":
           "mean immune interface", "mean Tcell NOT in interface": "mean immune NOT in interface",
           "T cell membrane ratio": "Immune membrane ratio"}
MACRO_NAMES = {DOUBLETS_TABLE: "Single Tumor Single T cell"}


def _read(folder, name):
    return pd.read_csv(Path(folder) / f"{name}.tsv", sep="\t")


def _read_reference(folder, name):
    return _read(folder, MACRO_NAMES.get(name, name)).rename(columns=RENAMED)


def test_cells_match_reference(result, reference):
    ref = _read_reference(reference, CELL_TABLE)
    own = _read(result.out_dir, CELL_TABLE)
    assert len(own) == len(ref)
    # The order of the cells within an event is not deterministic in Fiji, so the cells are matched on event and
    # area. The reference has no event column: it follows from the event table.
    events = _read_reference(reference, EVENT_TABLE)
    ref["Frame"] = np.repeat(events["Frame"], events["Immune cells"] + events["Tumor cells"]).to_numpy()
    m = ref.merge(own, on=["Frame", "Area"], suffixes=("_ref", "_own"), validate="one_to_one")
    assert len(m) == len(ref)
    for col in ("Immune intensity", "Tumor intensity", "Marker ratio", "Marker ratio over Area"):
        assert np.allclose(m[f"{col}_own"], m[f"{col}_ref"], atol=0.0015, rtol=0), col
    # Class here: 1 = T cell; in the macro an arbitrary cluster number, so compare the partition
    agree = (m["Class_own"] == m["Class_ref"]).mean()
    assert max(agree, 1 - agree) == 1.0


@pytest.mark.parametrize("name", [EVENT_TABLE, DOUBLETS_TABLE])
def test_event_tables_match_reference(result, reference, name):
    ref = _read_reference(reference, name)
    own = _read(result.out_dir, name)
    assert list(own.columns[: len(ref.columns)]) == list(ref.columns)  # the port adds columns at the end
    assert len(own) == len(ref)
    for col in ref.columns:
        assert np.allclose(own[col], ref[col], atol=0.0015, rtol=0, equal_nan=True), col


def test_membrane_distributions_match_reference(result, reference):
    """Not cell by cell: the port measures one membrane per cell (the macro one per ring of touching cells, I8) and
    a bootstrapped control (I10), so the distributions are compared."""
    ref = _read(reference, MEMBRANE_TABLE)
    own = _read(result.out_dir, MEMBRANE_TABLE)
    of_interest = _config()[0].channels.cell_of_interest
    wanted = "Immune" if of_interest == CellOfInterest.IMMUNE_CELL else "Tumor"
    assert len(own) == (result.cells["Cell type"] == wanted).sum()  # one row per cell of interest
    ref_if, own_if = ref["interface Mean"].notna(), own["interface Mean"].notna()
    assert abs(own_if.sum() - ref_if.sum()) / ref_if.sum() < 0.1
    for col, tol in [("membrane length", 0.05), ("non-interface Mean", 0.05), ("interface Mean", 0.05),
                     ("Ratio interface/non-interface", 0.05)]:
        assert own[col].median() == pytest.approx(ref[col].median(), rel=tol), col
    for q in (0.25, 0.5, 0.75):
        assert own["stdDevs above mean"].quantile(q) == pytest.approx(ref["stdDevs above mean"].quantile(q), abs=0.5)
    control = _read(result.out_dir, "Control windows stdDevs above mean").iloc[:, 0].dropna()
    for q in (0.25, 0.5, 0.75):
        assert control.quantile(q) == pytest.approx(ref["control stdDevs above mean"].quantile(q), abs=0.25)


def test_average_profile_matches_reference_cells_with_interface(result, reference):
    """The macro averages all rings (I11); recomputed from its profiles for the rings with an interface only."""
    import tifffile

    profiles = tifffile.imread(reference / "membrane_profiles_normalized.tif").astype(float)
    has = _read(reference, MEMBRANE_TABLE)["interface Mean"].notna().to_numpy()
    own_avg = _read(result.out_dir, AVERAGE_PROFILE)["profile"].to_numpy()
    n = min(len(own_avg), profiles.shape[1])
    ref_avg = np.nanmean(profiles[has][:, :n], axis=0)
    assert np.corrcoef(own_avg[:n], ref_avg)[0, 1] > 0.99
    assert own_avg[20] == pytest.approx(ref_avg[20], rel=0.1)
