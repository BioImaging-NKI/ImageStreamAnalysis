"""The user manual: the generated settings reference is up to date, and the site builds without warnings."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools" / "manual"))


def test_settings_reference_is_up_to_date():
    import gen_reference

    page = gen_reference.OUT.read_text(encoding="utf-8")
    assert page == gen_reference.render(), "run: python tools/manual/gen_reference.py"


def test_every_setting_is_described():
    import gen_reference

    source = gen_reference.CONFIG.read_text(encoding="utf-8")
    for _, cls, _ in gen_reference.SECTIONS:
        for name, _, comment in gen_reference.fields(source, cls):
            assert comment, f"{cls}.{name} has no comment in config.py (the manual's settings reference)"


@pytest.mark.skipif(shutil.which("mkdocs") is None and not (Path(sys.executable).parent / "mkdocs.exe").exists(),
                    reason="mkdocs not installed")
def test_manual_builds_strict(tmp_path):
    result = subprocess.run([sys.executable, "-m", "mkdocs", "build", "--strict", "--site-dir", str(tmp_path / "site")],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-3000:]
