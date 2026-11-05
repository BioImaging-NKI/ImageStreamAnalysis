"""Shortcuts to the import page with the logo as their icon (`imagestream-cli shortcut`; Windows).

The program file itself (`imagestream.exe` in the environment's Scripts folder) is a generic launcher that the
installer (uv, pip) writes again at every install, so it cannot keep an icon of its own. A shortcut can: it points
at that launcher and takes its icon from a copy of templates/logo.ico in the program's home folder.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

from .launcher import STATE_DIR

NAME = "ImageStream Analysis"


def launcher() -> Path:
    """imagestream.exe of this environment (next to python.exe, as the installer puts it)."""
    folder = Path(sys.executable).parent
    for candidate in (folder / "imagestream.exe", folder / "Scripts" / "imagestream.exe"):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"imagestream.exe not found next to {sys.executable}: install the package first")


def icon() -> Path:
    """The logo as an .ico in the home folder of the program (it outlives a reinstall of the package)."""
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    target = STATE_DIR / "logo.ico"
    with resources.as_file(resources.files("imagestream_analysis").joinpath("templates", "logo.ico")) as source:
        shutil.copyfile(source, target)
    return target


def _ps(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def make_shortcuts(desktop: bool = True, start_menu: bool = True) -> list[Path]:
    """Write the shortcuts; returns their paths."""
    if sys.platform != "win32":
        raise OSError("Shortcuts are made on Windows only")
    exe, ico = launcher(), icon()
    places = [p for p, wanted in (("Desktop", desktop), ("Programs", start_menu)) if wanted]
    script = "$s = New-Object -ComObject WScript.Shell\n"
    for place in places:
        script += (f"$p = Join-Path ([Environment]::GetFolderPath('{place}')) {_ps(NAME + '.lnk')}\n"
                   f"$l = $s.CreateShortcut($p); $l.TargetPath = {_ps(str(exe))}; "
                   f"$l.WorkingDirectory = {_ps(str(Path.home()))}; $l.IconLocation = {_ps(str(ico) + ',0')}; "
                   f"$l.Description = 'Analyse ImageStream data: the import and run page'; $l.Save(); $p\n")
    out = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", script],
                         capture_output=True, text=True, check=True)
    return [Path(line) for line in out.stdout.splitlines() if line.strip()]
