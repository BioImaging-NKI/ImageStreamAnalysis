import json
import threading
import tomllib
import urllib.request
from http.server import ThreadingHTTPServer

import numpy as np
import pandas as pd
import pytest

from conftest import IDEAS_IDS, write_ideas_zip
from imagestream_analysis import gui
from imagestream_analysis.config import ChannelSettings, Pool
from imagestream_analysis.pipeline import CELL_TABLE, EVENT_TABLE
from imagestream_analysis.project import OVERVIEW_FILE, DatasetEntry, Project, ZipEntry, run_project

OLD_ROLES = {"brightfield": 1, "marker_of_interest": 3, "t_cell": 6, "tumor": 9}  # projects before 2026-10-31
NAMES = ["BF", "MOI's name", "T cells", "Tumor"]


def make_project(ideas_zip, folder) -> Project:
    p = Project(folder)
    p.settings.segmentation.method = "labelmap"  # no Cellpose in the tests
    p.settings.segmentation.labelmap_file = ideas_zip["labelmap"]
    p.settings.membrane.thickness = 2
    p.settings.membrane.analysis = True
    ch = p.settings.channels  # the panel, by name
    ch.brightfield, ch.marker_of_interest, ch.immune, ch.tumor = "BF", ["MOI's name"], ["T cells"], "Tumor"
    p.zips.append(ZipEntry(ideas_zip["zip"], 0.5, [1, 3, 6, 9], list(NAMES), {},
                           [DatasetEntry("Sample A", "Experiment 1")]))
    return p


def test_panel_by_name_gives_the_numbers_of_each_zip(ideas_zip):
    p = make_project(ideas_zip, ideas_zip["tmp"] / "project")
    # a second export of the same panel, with its channels numbered differently
    p.zips.append(ZipEntry(ideas_zip["zip"], 0.5, [2, 4, 5, 8], ["Tumor", "BF", "MOI's name", "T cells"], {},
                           [DatasetEntry("Sample A", "Experiment 2")]))
    first, second = (p.channels_of(z) for z in p.zips)
    assert (first.brightfield, first.marker_of_interest, first.immune, first.tumor) == (1, [3], [6], 9)
    assert (second.brightfield, second.marker_of_interest, second.immune, second.tumor) == (4, [5], [8], 2)


def test_old_projects_with_roles_per_zip(ideas_zip):
    p = make_project(ideas_zip, ideas_zip["tmp"] / "project")
    p.settings.channels = ChannelSettings()
    p.zips[0].roles = dict(OLD_ROLES)
    text = p.to_toml()
    assert "t_cell = 6" in text
    loaded = Project.from_dict(tomllib.loads(text), p.folder)  # all zips agree: a panel by name
    ch = loaded.settings.channels
    assert (ch.brightfield, ch.marker_of_interest, ch.immune, ch.tumor) == ("BF", ["MOI's name"], ["T cells"],
                                                                           "Tumor")
    assert loaded.zips[0].roles == {} and loaded.problems() == []
    # an unnamed role channel: the roles stay as numbers, and still work
    data = tomllib.loads(text)
    data["zip"][0]["channel_names"][0] = ""
    loaded = Project.from_dict(data, p.folder)
    assert loaded.zips[0].roles == {"brightfield": 1, "marker_of_interest": [3], "immune": [6], "tumor": 9}
    assert loaded.channels_of(loaded.zips[0]).immune == [6]


def test_project_round_trip(ideas_zip):
    p = make_project(ideas_zip, ideas_zip["tmp"] / "project")
    path = p.save()
    loaded = Project.load(path.parent)
    assert loaded.zips == p.zips  # including a name with a quote
    assert loaded.settings == p.settings


def test_problems(ideas_zip):
    p = make_project(ideas_zip, ideas_zip["tmp"] / "project")
    assert p.problems() == []
    p.zips[0].pixel_size = None
    p.settings.channels.tumor = ""
    p.settings.channels.functional = ["CD-missing"]
    msgs = p.problems()
    assert any("pixel size" in m for m in msgs) and any("no channel for the tumor" in m for m in msgs)
    p.settings.channels.tumor = "Tumor"
    assert any("'CD-missing'" in m for m in p.problems())
    p.settings.channels.immune = []
    assert "No immune-cell marker" in p.problems()
    p.zips[0].datasets[0].include = False
    assert "No datasets selected" in p.problems()


def test_problems_list_a_membrane_pair_that_is_no_selection(ideas_zip):
    """A pair typed by hand (old projects) that matches no subset: listed, with the likely meaning."""
    p = make_project(ideas_zip, ideas_zip["tmp"] / "project")
    p.settings.channels.functional = ["T cells"]
    p.settings.membrane.cells = ["Immune T cells"]
    assert any("did you mean 'T cells+'" in m and "'Immune T cells+'" in m for m in p.problems())
    p.settings.membrane.cells = ["Immune T cells+"]
    assert p.problems() == []


def test_run_project_writes_results_and_overview(ideas_zip):
    folder = ideas_zip["tmp"] / "project"
    make_project(ideas_zip, folder).save()
    table = run_project(folder)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["experiment"] == "Experiment 1" and row["dataset"] == "Sample A"
    assert row["cells"] == 20 and row["immune cells"] == 10 and row["single tumor + single immune cell"] == 10
    assert row["with interface (MOI's name)"] == 10
    assert row["MOI on (MOI's name)"] == "Tumor" and row["membrane pair (MOI's name)"] == "Tumor with Immune"
    assert row["median Tumor membrane ratio, doublets (MOI's name)"] > 1
    out = folder / "Experiment 1" / "Sample A"
    cells = pd.read_csv(out / f"{CELL_TABLE}.tsv", sep="\t")
    assert sorted(set(cells["Object ID"])) == sorted(IDEAS_IDS)
    events = pd.read_csv(out / f"{EVENT_TABLE}.tsv", sep="\t")
    assert events["Object ID"].tolist() == sorted(IDEAS_IDS)
    assert (folder / OVERVIEW_FILE).exists()


def test_two_markers_of_interest_in_a_project(ideas_zip):
    """A second marker of interest (on the immune cells): its own pair from the data, its enrichment in the overview,
    and its rows in the project's membrane pairs.tsv."""
    from conftest import make_events

    data, _ = make_events(second_marker=True)
    path = write_ideas_zip(ideas_zip["tmp"] / "two.zip", data[:, [0, 1, 3, 2, 4]], IDEAS_IDS, [1, 3, 6, 9, 11])
    folder = ideas_zip["tmp"] / "project"
    p = make_project(ideas_zip, folder)
    p.zips[0] = ZipEntry(path, 0.5, [1, 3, 6, 9, 11], [*NAMES, "Second"], {},
                         [DatasetEntry("Sample A", "Experiment 1")])
    p.settings.channels.marker_of_interest = ["MOI's name", "Second"]
    p.save()
    assert Project.load(folder).settings.channels.marker_of_interest == ["MOI's name", "Second"]
    row = run_project(folder).iloc[0]
    assert row["MOI on (MOI's name)"] == "Tumor" and row["membrane pair (MOI's name)"] == "Tumor with Immune"
    assert row["MOI on (Second)"] == "Immune" and row["membrane pair (Second)"] == "Immune with Tumor"
    assert row["median stdDevs above mean (Second)"] > 1 and row["median Immune membrane ratio, doublets (Second)"] > 1
    pairs = pd.read_csv(folder / "membrane pairs.tsv", sep="\t")
    assert list(pairs.columns[:4]) == ["experiment", "dataset", "marker", "marker on"]
    assert pairs.groupby("marker")["marker on"].first().to_dict() == {"MOI's name": "Tumor", "Second": "Immune"}
    assert (folder / "Experiment 1" / "Sample A" / "Membrane intensity analysis - Second.tsv").exists()


def test_older_results_in_the_overview(ideas_zip):
    """A dataset folder written before 2026-11-03 (Membrane pairs.tsv without marker columns) is still read."""
    from imagestream_analysis.project import write_overview

    folder = ideas_zip["tmp"] / "project"
    make_project(ideas_zip, folder).save()
    run_project(folder)
    f = folder / "Experiment 1" / "Sample A" / "Membrane pairs.tsv"
    pd.read_csv(f, sep="\t").drop(columns=["marker", "marker on"]).to_csv(f, sep="\t", index=False)
    assert len(write_overview(folder)) == 1
    pairs = pd.read_csv(folder / "membrane pairs.tsv", sep="\t", keep_default_na=False)
    assert (pairs["marker"] == "").all() and "membrane of" in pairs


def test_markers_of_interest_from_the_page(ideas_zip):
    """The page sends a list; an empty field (None or "") means no marker of interest, not the default channel."""
    data = gui.project_to_json(make_project(ideas_zip, ideas_zip["tmp"] / "p"))
    assert data["settings"]["channels"]["marker_of_interest"] == ["MOI's name"]
    for value, expected in ((["MOI's name", "T cells"], ["MOI's name", "T cells"]), ([], []), (None, []), ("", [])):
        data["settings"]["channels"]["marker_of_interest"] = value
        assert gui.project_from_json(data).settings.channels.marker_of_interest == expected


def test_a_failing_dataset_does_not_stop_the_others(ideas_zip):
    folder = ideas_zip["tmp"] / "project"
    p = make_project(ideas_zip, folder)
    p.zips[0].datasets.append(DatasetEntry("missing", "Experiment 2"))
    p.save()
    table = run_project(folder)
    assert len(table) == 2
    failed = table[table["dataset"] == "missing"].iloc[0]
    assert "no dataset" in failed["error"]


@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "STATE_DIR", tmp_path / "state")
    monkeypatch.setattr(gui, "PREFS", tmp_path / "state" / "prefs.json")
    token = "test-token"
    srv = ThreadingHTTPServer(("127.0.0.1", 0), gui.make_handler(token, gui.Job(), [0.0]))
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}", token
    srv.shutdown()


def call(server, api, data=None, token=None):
    url, tok = server
    req = urllib.request.Request(f"{url}/api/{api}", data=json.dumps(data or {}).encode(), method="POST",
                                 headers={"X-Token": token if token is not None else tok})
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())


def test_page_and_token(server):
    url, token = server
    with urllib.request.urlopen(url + "/") as r:
        html = r.read().decode()
    assert token in html and "ImageStream Analysis" in html
    with pytest.raises(urllib.error.HTTPError) as err:
        call(server, "state", token="wrong")
    assert err.value.code == 403


def test_ping_keeps_alive_and_quit_stops(tmp_path, monkeypatch):
    """The page's ping counts as use (the server stops IDLE_TIMEOUT after the last one); Quit stops the server."""
    monkeypatch.setattr(gui, "STATE_DIR", tmp_path / "state")
    last_seen = [0.0]
    srv = ThreadingHTTPServer(("127.0.0.1", 0), gui.make_handler("t", gui.Job(), last_seen))
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{srv.server_address[1]}"
    with urllib.request.urlopen(url + "/ping") as r:
        assert r.read() == b"imagestream-analysis"
    assert last_seen[0] > 0
    with urllib.request.urlopen(url + "/") as r:
        assert 'id="quit"' in r.read().decode()
    assert call((url, "t"), "quit")["ok"]
    thread.join(timeout=5)
    assert not thread.is_alive()
    srv.server_close()


def test_api_scan_save_and_overview(server, ideas_zip):
    state = call(server, "state")["result"]
    assert "brightfield" in state["roles"] and state["defaults"]["segmentation"]["model"] == "cyto3"
    scanned = call(server, "scan", {"path": str(ideas_zip["zip"])})["result"]
    assert scanned["datasets"][0]["objects"] == len(IDEAS_IDS)
    assert sorted(int(c) for c in scanned["thumbnails"]) == [1, 3, 6, 9]
    bad = call(server, "scan", {"path": str(ideas_zip["tmp"] / "nothing.zip")})
    assert not bad["ok"] and "zip" in bad["error"]

    folder = ideas_zip["tmp"] / "from page"
    project = gui.project_to_json(make_project(ideas_zip, folder))
    saved = call(server, "save", {"project": project})["result"]
    assert saved["problems"] == [] and (folder / "project.toml").exists()
    assert Project.load(folder).settings.channels.immune == ["T cells"]
    run_project(folder)
    o = call(server, "overview", {"folder": str(folder)})["result"]
    assert o["rows"][0]["cells"] == 20 and "experiment" in o["columns"]
    unknown = call(server, "run", {"project": project, "datasets": ["Nothing"]})  # run one dataset: checked first
    assert not unknown["ok"] and "No included dataset named Nothing" in unknown["error"]


def test_apply_from_a_dashboard(ideas_zip, monkeypatch):
    """Apply: moved thresholds (None = automatic again) and the membrane pair into the project, then a run."""
    folder = ideas_zip["tmp"] / "applied"
    project = make_project(ideas_zip, folder)
    project.settings.classification.thresholds = {"Old": 1.0, "Kept": 2.0}
    project.save()
    monkeypatch.setitem(gui.SERVED, "root", folder)

    class FakeJob:
        running, started = False, None

        def start(self, path):
            self.started = path

    job = FakeJob()
    assert gui.apply_from_dashboard({"dry_run": True}, job) == {"datasets": 1, "running": False}
    result = gui.apply_from_dashboard({"thresholds": {"A": 3.5, "Old": None},
                                       "membrane": {"cells": ["Tumor"], "partners": ["A+"]}}, job)
    assert result["started"] and job.started == folder / "project.toml"
    s = Project.load(folder).settings
    assert s.classification.thresholds == {"Kept": 2.0, "A": 3.5}
    assert s.membrane.cells == ["Tumor"] and s.membrane.partners == ["A+"]
    monkeypatch.setitem(gui.SERVED, "root", ideas_zip["tmp"])  # no project there
    with pytest.raises(ValueError, match="does not belong"):
        gui.apply_from_dashboard({"dry_run": True}, job)


def test_object_ids_follow_the_frames():
    from imagestream_analysis.pipeline import _with_object_ids

    t = pd.DataFrame({"Frame": [2, 1, 3]})
    assert _with_object_ids(t, np.array([7, 9, 11]))["Object ID"].tolist() == [9, 7, 11]
    assert "Object ID" not in _with_object_ids(t, None)


@pytest.mark.parametrize("value", [None, ""])
def test_empty_max_events_from_the_page_means_all_events(ideas_zip, value):
    data = gui.project_to_json(make_project(ideas_zip, ideas_zip["tmp"] / "p"))
    assert data["settings"]["input"]["max_events"] is None  # the default: all events
    data["settings"]["input"]["max_events"] = 500
    assert gui.project_from_json(data).settings.input.max_events == 500
    data["settings"]["input"]["max_events"] = value  # the field cleared on the page
    p = gui.project_from_json(data)
    assert p.settings.input.max_events is None
    p.save()
    assert Project.load(p.folder).settings.input.max_events is None
    assert all(s.input.max_events is None for _, _, s in p.samples())


@pytest.mark.parametrize("pool", ["experiment", "dataset"])
def test_datasets_of_an_experiment_are_classified_together(ideas_zip, pool):
    folder = ideas_zip["tmp"] / "project"
    p = make_project(ideas_zip, folder)
    p.settings.classification.pool = Pool(pool)
    p.settings.channels.functional = ["MOI's name"]
    second = write_ideas_zip(ideas_zip["tmp"] / "second.zip", ideas_zip["data"][:, [0, 1, 2, 3]], IDEAS_IDS,
                             [1, 3, 6, 9], dataset="Sample B")
    p.zips.append(ZipEntry(second, 0.5, [1, 3, 6, 9], list(NAMES), {}, [DatasetEntry("Sample B", "Experiment 1")]))
    p.save()
    table = run_project(folder)
    assert len(table) == 2 and (table["cells"] == 20).all()
    a, b = (folder / "Experiment 1" / name for name in ("Sample A", "Sample B"))
    ta, tb = (pd.read_csv(f / "Marker thresholds.tsv", sep="\t") for f in (a, b))
    info = json.loads((a / "run.json").read_text(encoding="utf-8"))
    if pool == "experiment":
        assert ta.equals(tb) and ta["of"].iloc[0] == 20  # one set of thresholds for the immune cells of both
        assert info["classified with"] == "1 other dataset(s)"
    else:
        assert ta["of"].iloc[0] == 10 and info["classified with"] == "this dataset only"
    inter = pd.read_csv(folder / "interactions.tsv", sep="	")
    total = inter[(inter["dataset"] == "(all)") & (inter["cell type"] == "Tumor")
                  & inter["neighbour type"].str.startswith("Immune")]  # immune subtypes: one row each
    assert total["contacts"].sum() == 20 and (total["cells"] == 20).all()  # both datasets of the experiment
    assert (table["interactions"] == 10).all()
    log_text = (a / "log.txt").read_text(encoding="utf-8")
    assert "Marker normalization" in log_text and "Done in" in log_text  # one log over both stages


def test_one_immune_marker_per_zip_is_enough():
    """The immune markers are summed: a zip needs one of them, not all (another zip may name its marker otherwise)."""
    ch = ChannelSettings(brightfield="BF", marker_of_interest=[], tumor="Tumor", immune=["CD3", "T cells"],
                         functional=[])
    numbers, names = [1, 2, 3], ["BF", "Tumor", "T cells"]
    assert ch.resolved(numbers, names).immune == [3]
    with pytest.raises(ValueError, match="immune marker 'CD3'"):
        ch.resolved(numbers, ["BF", "Tumor", "other"])


@pytest.mark.parametrize("pool", ["experiment", "dataset"])
def test_run_one_dataset_again(ideas_zip, pool):
    """Run only one dataset: the other keeps its results (with pools per experiment its experiment runs along)."""
    folder = ideas_zip["tmp"] / "project"
    p = make_project(ideas_zip, folder)
    p.settings.classification.pool = Pool(pool)
    second = write_ideas_zip(ideas_zip["tmp"] / "second.zip", ideas_zip["data"][:, [0, 1, 2, 3]], IDEAS_IDS,
                             [1, 3, 6, 9], dataset="Sample B")
    p.zips.append(ZipEntry(second, 0.5, [1, 3, 6, 9], list(NAMES), {}, [DatasetEntry("Sample B", "Experiment 1")]))
    p.save()
    run_project(folder)
    a, b = (folder / "Experiment 1" / name / "run.json" for name in ("Sample A", "Sample B"))
    before = a.stat().st_mtime_ns, b.stat().st_mtime_ns
    table = run_project(folder, ["Sample B"])
    assert len(table) == 2  # the overview still has both
    assert b.stat().st_mtime_ns != before[1]
    assert (a.stat().st_mtime_ns != before[0]) == (pool == "experiment")  # classified together: runs along
    with pytest.raises(ValueError, match="No included dataset named Nothing"):
        run_project(folder, ["Nothing"])
