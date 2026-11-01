"""The import and run page: a small local web app (`imagestream-cli gui`, or `imagestream`).

Add IDEAS exports, name the channels and give them their role, set the pixel size and the analysis settings, run
all datasets and see the overview of the results. The project is saved as `project.toml` in its output folder and
run as a subprocess (`imagestream-cli run-project`), whose log is streamed to the page.

The server only listens on 127.0.0.1 and every API call must carry a random per-session token, so other web pages
in the browser cannot use it. The page pings it every minute; it stops IDLE_TIMEOUT after the last ping or request
(so 15 minutes after the page was closed; never during a run), or at once with the page's Quit button.
"""

from __future__ import annotations

import json
import logging
import mimetypes
import os
import secrets
import subprocess
import sys
import threading
import time
import urllib.request
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__, ideas
from .config import THRESHOLD_METHODS, CellOfInterest, ChannelSettings, MarkerLocation, Pool, Settings
from .dashboard import page_template
from .project import INTERACTIONS_FILE, Project, overview

log = logging.getLogger(__name__)
STATE_DIR = Path(os.environ.get("ISA_HOME") or Path.home() / ".imagestream_analysis")  # ISA_HOME: a test instance
INSTANCE = STATE_DIR / "gui_instance.json"  # port of a running page server, so that a second start reuses it
PREFS = STATE_DIR / "gui_prefs.json"  # last project folder, light / dark
IDLE_TIMEOUT = 15 * 60  # s after the last request; the page pings every PING s
PING = 60  # s
NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)


def load_prefs() -> dict[str, Any]:
    try:
        prefs: dict[str, Any] = json.loads(PREFS.read_text(encoding="utf-8"))
        return prefs
    except (OSError, ValueError):
        return {}


def save_prefs(changes: dict[str, Any]) -> dict[str, Any]:
    prefs = load_prefs() | changes
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    PREFS.write_text(json.dumps(prefs, indent=2), encoding="utf-8")
    return prefs


# ---- the project as the page sees it (JSON) <-> Project


def project_to_json(p: Project) -> dict[str, Any]:
    return {"folder": str(p.folder), "settings": p.settings.to_dict(),
            "zips": [{"path": str(z.path), "pixel_size": z.pixel_size, "channels": z.channels,
                      "channel_names": z.channel_names, "roles": z.roles,
                      "datasets": [{"name": d.name, "experiment": d.experiment, "include": d.include}
                                   for d in z.datasets]} for z in p.zips]}


def project_from_json(data: dict[str, Any]) -> Project:
    folder = Path(data.get("folder") or "")
    if not str(data.get("folder") or "").strip():
        raise ValueError("Choose an output folder for the project")
    settings = {k: {kk: vv for kk, vv in v.items() if vv is not None and vv != ""}
                for k, v in (data.get("settings") or {}).items()}
    # a role without a channel stays empty (not the default channel number of a settings file)
    channels = (data.get("settings") or {}).get("channels") or {}
    for role in ("brightfield", "marker_of_interest", "tumor"):
        if role in channels and channels[role] in (None, ""):
            settings.setdefault("channels", {})[role] = ""
    zips = [z | {"dataset": z.get("datasets", [])} for z in data.get("zips", [])]
    return Project.from_dict(settings | {"zip": zips}, folder)


# ---- requests


def scan_zip(path: str) -> dict[str, Any]:
    """The datasets and channels of an IDEAS zip, with example images of the first dataset."""
    p = Path(path)
    if p.suffix.lower() != ".zip" or not p.is_file():
        raise ValueError(f"Not a zip file: {path}")
    datasets = ideas.scan(p)
    if not datasets:
        raise ValueError(f"{p.name} holds no IDEAS export (folders with <object>_Ch<N>.ome.tif files)")
    first = next(iter(datasets))
    return {"path": str(p), "name": p.name,
            "datasets": [{"name": d.name, "objects": len(d.files), "complete": len(d.objects()), "channels": d.channels}
                         for d in datasets.values()],
            "thumbnails": {str(c): url for c, url in ideas.thumbnails(p, first).items()}}


def browse(kind: str, initial: str = "") -> list[str]:
    """A native file or folder dialog (tkinter), in its own thread so that the server keeps running."""
    result: list[str] = []

    def ask() -> None:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        start = initial if initial and Path(initial).exists() else str(Path.home())
        if kind == "zip":
            result.extend(filedialog.askopenfilenames(parent=root, title="IDEAS exports", initialdir=start,
                                                      filetypes=[("IDEAS export", "*.zip")]))
        elif kind == "project":
            f = filedialog.askopenfilename(parent=root, title="Open project", initialdir=start,
                                           filetypes=[("Project", "project.toml")])
            if f:
                result.append(f)
        else:
            f = filedialog.askdirectory(parent=root, title="Output folder of the project", initialdir=start)
            if f:
                result.append(f)
        root.destroy()

    t = threading.Thread(target=ask)
    t.start()
    t.join()
    return [str(Path(r)) for r in result]


class Job:
    """One run of `imagestream-cli run-project` at a time; its output is kept for the page."""

    def __init__(self) -> None:
        self.proc: subprocess.Popen[str] | None = None
        self.lines: list[str] = []
        self.progress = ""
        self.folder = ""
        self.lock = threading.Lock()

    @property
    def running(self) -> bool:
        return self.proc is not None and self.proc.poll() is None

    def start(self, project_file: Path) -> None:
        if self.running:
            raise RuntimeError("A run is still going")
        exe = Path(sys.executable)
        python = exe.with_name("python.exe") if exe.name.lower() == "pythonw.exe" else exe
        cmd = [str(python), "-m", "imagestream_analysis.cli", "run-project", str(project_file)]
        with self.lock:
            self.lines, self.progress, self.folder = [], "", str(project_file.parent)
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                     encoding="utf-8", errors="replace", creationflags=NO_WINDOW,
                                     env=os.environ | {"PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
        threading.Thread(target=self._read, daemon=True).start()

    def _read(self) -> None:
        assert self.proc is not None and self.proc.stdout is not None
        for line in self.proc.stdout:
            line = line.rstrip()
            with self.lock:
                if line.startswith("ISA_PROGRESS "):
                    self.progress = line[len("ISA_PROGRESS "):]
                else:
                    self.lines.append(line)

    def stop(self) -> None:
        if self.running and self.proc is not None:
            self.proc.terminate()

    def status(self, since: int) -> dict[str, Any]:
        with self.lock:
            return {"running": self.running, "lines": self.lines[since:], "n": len(self.lines),
                    "progress": self.progress, "folder": self.folder,
                    "returncode": None if self.proc is None else self.proc.poll()}


SERVED: dict[str, Path] = {}  # the project folder whose files (dashboards) the page may open
CONTENT_TYPES = {".webp": "image/webp", ".html": "text/html; charset=utf-8", ".tsv": "text/tab-separated-values"}


def served_file(url_path: str) -> Path | None:
    """A file under the served project folder, or None (outside it, or not there)."""
    root = SERVED.get("root")
    if root is None:
        return None
    from urllib.parse import unquote

    target = (root / unquote(url_path)).resolve()
    if not target.is_file() or not target.is_relative_to(root.resolve()):
        return None
    return target


def handle(api: str, data: dict[str, Any], job: Job) -> Any:
    if api == "state":
        prefs = load_prefs()
        project = None
        last = prefs.get("project")
        if last and Path(last, "project.toml").exists():
            try:
                project = project_to_json(Project.load(last))
            except Exception as e:  # noqa: BLE001 - a broken file must not stop the page
                log.warning("Could not open the last project %s: %s", last, e)
        return {"version": __version__, "theme": prefs.get("theme", ""), "project": project,
                "defaults": Settings().to_dict(), "roles": list(ChannelSettings().roles()),
                "marker_locations": [str(m) for m in MarkerLocation],
                "cells_of_interest": [str(c) for c in CellOfInterest],
                "threshold_methods": THRESHOLD_METHODS, "pools": [str(p) for p in Pool]}
    if api == "prefs":
        return save_prefs({k: v for k, v in data.items() if k in ("theme",)})
    if api == "browse":
        return browse(str(data.get("kind", "zip")), str(data.get("initial", "")))
    if api == "scan":
        return scan_zip(str(data["path"]))
    if api == "open_project":
        p = Project.load(str(data["path"]))
        save_prefs({"project": str(p.folder)})
        return project_to_json(p)
    if api == "save":
        p = project_from_json(data["project"])
        path = p.save()
        save_prefs({"project": str(p.folder)})
        return {"path": str(path), "problems": p.problems()}
    if api == "run":
        p = project_from_json(data["project"])
        problems = p.problems()
        if problems:
            return {"started": False, "problems": problems}
        path = p.save()
        save_prefs({"project": str(p.folder)})
        job.start(path)
        return {"started": True, "problems": []}
    if api == "status":
        return job.status(int(data.get("since", 0)))
    if api == "stop":
        job.stop()
        return {"stopped": True}
    if api == "overview":
        folder = Path(str(data["folder"]))
        if folder.is_dir():
            SERVED["root"] = folder
        table = overview(folder) if folder.is_dir() else None
        if table is None or table.empty:
            return {"columns": [], "rows": [], "interactions": []}
        contacts = folder / INTERACTIONS_FILE
        summed = pd.read_csv(contacts, sep="	") if contacts.exists() else pd.DataFrame()
        if len(summed):
            summed = summed[summed["dataset"] == "(all)"]
        return {"columns": list(table.columns),
                "rows": json.loads(table.to_json(orient="records", force_ascii=False)),
                "interactions": json.loads(summed.to_json(orient="records", force_ascii=False)) if len(summed) else []}
    if api == "open":
        target = Path(str(data["path"]))
        if not target.exists():
            raise ValueError(f"Not found: {target}")
        os.startfile(target)  # Windows only, as the rest of the page
        return {"opened": str(target)}
    raise ValueError(f"Unknown request: {api}")


def page(token: str) -> bytes:
    html = page_template("gui.html")
    return (html.replace("__ISA_TOKEN__", token).replace("__ISA_VERSION__", __version__)
            .replace("__ISA_PING__", str(PING * 1000)).encode("utf-8"))


def make_handler(token: str, job: Job, last_seen: list[float]) -> type[BaseHTTPRequestHandler]:
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002 - quiet console
            pass

        def _send(self, code: int, body: bytes, ctype: str) -> None:
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self) -> None:  # noqa: N802
            last_seen[0] = time.monotonic()
            if self.path in ("/", "/index.html"):
                self._send(200, page(token), "text/html; charset=utf-8")
            elif self.path.startswith(f"/files/{token}/"):
                target = served_file(self.path[len(f"/files/{token}/"):].split("?")[0])
                if target is None:
                    self._send(404, b"not found", "text/plain")
                    return
                ctype = CONTENT_TYPES.get(target.suffix.lower()) or mimetypes.guess_type(target.name)[0]
                self._send(200, target.read_bytes(), ctype or "application/octet-stream")
            elif self.path == "/ping":
                self._send(200, b"imagestream-analysis", "text/plain")
            else:
                self._send(404, b"not found", "text/plain")

        def do_POST(self) -> None:  # noqa: N802
            last_seen[0] = time.monotonic()
            if not self.path.startswith("/api/") or self.headers.get("X-Token") != token:
                self._send(403, b"forbidden", "text/plain")
                return
            if self.path == "/api/quit":  # the Quit button: stop a run, answer, then stop the server
                job.stop()
                self._send(200, b'{"ok": true, "result": null}', "application/json")
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            try:
                length = int(self.headers.get("Content-Length", 0))
                data = json.loads(self.rfile.read(length) or b"{}")
                body = {"ok": True, "result": handle(self.path[len("/api/"):], data, job)}
            except Exception as e:  # noqa: BLE001 - reported on the page
                log.warning("%s failed: %s", self.path, e)
                body = {"ok": False, "error": str(e)}
            self._send(200, json.dumps(body, default=str).encode("utf-8"), "application/json")

    return Handler


def _running_instance() -> str | None:
    try:
        url = json.loads(INSTANCE.read_text(encoding="utf-8"))["url"]
        with urllib.request.urlopen(url.rstrip("/") + "/ping", timeout=1) as r:
            return str(url) if r.read() == b"imagestream-analysis" else None
    except Exception:  # noqa: BLE001 - no (live) instance
        return None


def serve(port: int = 0, open_browser: bool = True) -> None:
    url = _running_instance()
    if url and port == 0:  # a page server is running: show it again instead of starting a second one
        if open_browser:
            webbrowser.open(url)
        print(f"Already running: {url}")
        return
    token = secrets.token_urlsafe(24)
    job = Job()
    last_seen = [time.monotonic()]
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(token, job, last_seen))
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    INSTANCE.write_text(json.dumps({"url": url}), encoding="utf-8")

    def watchdog() -> None:
        while True:
            time.sleep(30)
            if not job.running and time.monotonic() - last_seen[0] > IDLE_TIMEOUT:
                server.shutdown()
                return

    threading.Thread(target=watchdog, daemon=True).start()
    print(f"ImageStream analysis: {url} (stops {IDLE_TIMEOUT // 60} min after the page was closed, or with Quit)")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    finally:
        server.server_close()
        try:
            if json.loads(INSTANCE.read_text(encoding="utf-8")).get("url") == url:
                INSTANCE.unlink()
        except (OSError, ValueError):
            pass


def main_gui() -> int:
    serve()
    return 0
