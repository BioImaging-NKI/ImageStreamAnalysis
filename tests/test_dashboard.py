import base64
import json
import re
import shutil
import subprocess
from importlib import resources

import numpy as np
import pytest

from imagestream_analysis import gui
from imagestream_analysis.dashboard import FILES, _drop_collinear, _outline
from imagestream_analysis.pipeline import run
from test_pipeline import settings_for


def embedded_data(html: str) -> dict:
    m = re.search(r"const D = (\{.*?\});\n", html, re.S)
    assert m, "no data in the page"
    return json.loads(m.group(1))


def test_dashboard_is_written_with_data_and_image_sheets(synthetic):
    s = settings_for(synthetic)
    s.membrane.analysis = True
    result = run(s)
    page = result.out_dir / "dashboard.html"
    data = embedded_data(page.read_text(encoding="utf-8"))
    assert data["n_events"] == 10 and len(data["cells"]["type"]) == 20
    assert data["tile"] == [60, 60]
    assert {c["role"] for c in data["channels"]} == {"brightfield", "marker_of_interest", "tumor", "immune"}
    for ch in data["channels"]:
        assert (result.out_dir / FILES / f"c{ch['number']}_0.webp").exists()
    outline = np.array(data["cells"]["outline"][0]).reshape(-1, 2) / 2
    assert outline.min() >= 0 and outline.max() <= 60
    assert all(len(o) >= 8 for o in data["cells"]["outline"])
    assert sum(len(i) for i in data["events"]["interfaces"]) >= 10  # every event has an interface
    m = data["membranes"]  # every cell's membrane; the page computes the enrichment of any pair
    assert sorted(m["cell"]) == list(range(20)) and len(m["markers"]) == 1
    first = m["markers"][0]  # the pair of the run: the marker is on the tumor cells
    assert (first["name"], first["carried_by"], first["cells"], first["partners"]) == ("MOI", "Tumor", ["Tumor"],
                                                                                        ["Immune"])
    profile = np.frombuffer(base64.b64decode(first["profile"]), "<f4")
    track = np.frombuffer(base64.b64decode(m["track"]), "<i4")
    assert len(profile) == len(track) == sum(m["length"])
    assert set(np.unique(track)) <= set(range(-1, 20)) and (track >= 0).any()
    assert first["expression"]["check"]
    assert len(data["pairs"]) == 20 and data["markers"] == []  # 10 touching pairs, as rows of the cell table


def test_dashboard_data_with_two_markers_of_interest(synthetic_two_markers):
    """The tracks once, a profile per marker, each with the pair of the run (the type that carries it)."""
    s = settings_for(synthetic_two_markers)
    s.channels.marker_of_interest = [2, 5]
    result = run(s)
    data = embedded_data((result.out_dir / "dashboard.html").read_text(encoding="utf-8"))
    m = data["membranes"]
    assert [(k["name"], k["carried_by"], k["cells"], k["partners"]) for k in m["markers"]] == [
        ("MOI", "Tumor", ["Tumor"], ["Immune"]), ("Second", "Immune", ["Immune"], ["Tumor"])]
    track = np.frombuffer(base64.b64decode(m["track"]), "<i4")
    for k in m["markers"]:
        assert len(np.frombuffer(base64.b64decode(k["profile"]), "<f4")) == len(track) == sum(m["length"])
    # the macro's per-event measure, per marker, on the side of the type that carries it
    assert [(e["name"], e["carried_by"]) for e in data["events"]["moi"]] == [("MOI", "Tumor"), ("Second", "Immune")]
    assert len(data["events"]["moi"][1]["t_ratio"]) == 10
    assert "iface-card" not in (result.out_dir / "dashboard.html").read_text(encoding="utf-8")


def test_outline_of_a_square_is_its_four_corners_on_pixel_edges():
    mask = np.zeros((10, 10), bool)
    mask[2:5, 3:7] = True  # rows 2-4, columns 3-6: edges x 3..7, y 2..5
    pts = np.array(_outline(mask)).reshape(-1, 2) / 2
    assert pts[:, 0].min() == pytest.approx(2.5, abs=0.5) and pts[:, 0].max() == pytest.approx(6.5, abs=0.5)
    assert pts[:, 1].min() == pytest.approx(1.5, abs=0.5) and pts[:, 1].max() == pytest.approx(4.5, abs=0.5)
    assert len(pts) <= 8  # straight runs are reduced to their ends


def test_drop_collinear():
    square = np.array([[0, 0], [1, 0], [2, 0], [2, 1], [2, 2], [2, 2], [1, 2], [0, 2], [0, 1]])
    assert _drop_collinear(square).tolist() == [[0, 0], [2, 0], [2, 2], [0, 2]]


def test_served_files_stay_inside_the_project(tmp_path, monkeypatch):
    root = tmp_path / "project"
    (root / "exp" / "ds").mkdir(parents=True)
    (root / "exp" / "ds" / "dashboard.html").write_text("x")
    (tmp_path / "secret.txt").write_text("no")
    monkeypatch.setitem(gui.SERVED, "root", root)
    assert gui.served_file("exp/ds/dashboard.html") == (root / "exp" / "ds" / "dashboard.html").resolve()
    assert gui.served_file("exp%2Fds%2Fdashboard.html") is not None
    assert gui.served_file("../secret.txt") is None
    assert gui.served_file("..%2Fsecret.txt") is None
    assert gui.served_file("exp/ds/missing.html") is None


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
@pytest.mark.parametrize("name", ["gui.html", "dashboard.html"])
def test_page_scripts_parse(name, tmp_path):
    html = resources.files("imagestream_analysis").joinpath("templates", name).read_text(encoding="utf-8")
    js = "\n".join(re.findall(r"<script>(.*?)</script>", html, re.S))
    src = tmp_path / "page.js"
    src.write_text(js, encoding="utf-8")
    proc = subprocess.run(["node", "--check", str(src)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
