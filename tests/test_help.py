"""The context help (templates/help.json, from the manual): up to date, and every ? button points at a section of it."""

import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "imagestream_analysis" / "templates"


def _help() -> dict:
    return json.loads((TEMPLATES / "help.json").read_text(encoding="utf-8"))


def test_help_is_up_to_date_with_the_manual():
    pytest.importorskip("pymdownx", reason="needs the docs dependencies (uv sync --group docs)")
    spec = importlib.util.spec_from_file_location("make_help", ROOT / "tools" / "manual" / "make_help.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert (TEMPLATES / "help.json").read_text(encoding="utf-8") == mod.text(), \
        "the manual changed: run python tools/manual/make_help.py"


def test_every_help_target_exists():
    """data-help, data-help-btn, data-help-scope, data-help-home and the targets in the scripts ('guide/x#anchor')."""
    sections = {p["path"]: set(re.findall(r'<section data-a="([^"]*)"', p["html"])) for p in _help()["pages"]}
    targets = set()
    for f in TEMPLATES.glob("*.html"):
        targets |= {(f.name, t) for t in re.findall(r'"((?:guide|analysis|reference|help)/[a-z-]+(?:#[a-z0-9-]*)?)"',
                                                    f.read_text(encoding="utf-8"))}
    assert len(targets) > 20
    bad = [(name, t) for name, t in sorted(targets)
           if t.split("#")[0] not in sections or ("#" in t and t.split("#")[1] not in sections[t.split("#")[0]])]
    assert not bad, f"help targets that are not in the manual (renamed headings?): {bad}"


def test_help_is_embedded_in_the_pages():
    from imagestream_analysis import gui
    from imagestream_analysis.dashboard import page_template

    for page in (gui.page("t").decode("utf-8"), page_template("dashboard.html")):
        assert "__ISA_HELP" not in page
        assert "#helpPanel {" in page and "const HELP" in page and 'id="helpBtn"' in page
        raw = re.search(r'<script id="isa-help" type="application/json">(.*?)</script>', page, re.DOTALL).group(1)
        assert json.loads(raw)["pages"]
