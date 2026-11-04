"""The import and run page: a small local web app (`imagestream-cli gui`, or `imagestream`).

Add IDEAS exports, name the channels and give them their role, set the pixel size and the analysis settings, run
all datasets and see the overview of the results. The project is saved as `project.toml` in its output folder and
run as a subprocess (`imagestream-cli run-project`), whose log is streamed to the page.

The server only listens on 127.0.0.1 and every API call must carry a random per-session token, so other web pages
in the browser cannot use it. The page pings it every minute; it stops IDLE_TIMEOUT (launcher.py) after the last ping
or request (so 15 minutes after the page was closed; never during a run), or at once with the page's Quit button.
"""

from __future__ import annotations

import copy
import json
import logging
import mimetypes
import os
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__, ideas
from .comparison import compare, write_comparison
from .config import THRESHOLD_METHODS, BandCentre, ChannelSettings, MarkerLocation, Pool, Settings
from .dashboard import page_template
from .launcher import STATE_DIR
from .project import INTERACTIONS_FILE, PROJECT_FILE, Project, in_project, overview, result_folders, select_samples

log = logging.getLogger(__name__)
PREFS = STATE_DIR / "gui_prefs.json"  # last project folder, light / dark
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
    """The project for the page. A membrane pair for every marker (membrane.cells / partners, older projects) is shown
    and saved as an entry of membrane.pairs per marker of interest: the page has a pair per marker only."""
    p = copy.deepcopy(p)
    p.pairs_per_marker()
    return {"folder": str(p.folder), "settings": p.settings.to_dict(),
            "zips": [{"path": str(z.path), "pixel_size": z.pixel_size, "channels": z.channels,
                      "channel_names": z.channel_names, "roles": z.roles,
                      "datasets": [{"name": d.name, "experiment": d.experiment, "include": d.include,
                                    "condition": d.condition, "replicate": d.replicate}
                                   for d in z.datasets]} for z in p.zips],
            "reference": p.reference}


def project_from_json(data: dict[str, Any]) -> Project:
    folder = Path(data.get("folder") or "")
    if not str(data.get("folder") or "").strip():
        raise ValueError("Choose an output folder for the project")
    settings = {k: {kk: vv for kk, vv in v.items() if vv is not None and vv != ""}
                for k, v in (data.get("settings") or {}).items()}
    # a role without a channel stays empty (not the default channel number of a settings file); the markers of
    # interest are a list (an empty one is kept above), "" from an older page means none
    channels = (data.get("settings") or {}).get("channels") or {}
    for role in ("brightfield", "marker_of_interest", "tumor"):
        if role in channels and channels[role] in (None, ""):
            settings.setdefault("channels", {})[role] = ""
    zips = [z | {"dataset": z.get("datasets", [])} for z in data.get("zips", [])]
    return Project.from_dict(settings | {"zip": zips, "comparison": {"reference": data.get("reference") or ""}},
                             folder)


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


def run_folders(project_file: Path, only: list[str] | None = None) -> dict[str, tuple[str, str]]:
    """The result folders a run of the project writes (relative, posix) -> (experiment, dataset)."""
    p = Project.load(project_file)
    folders = result_folders(p)
    if not only:
        return folders
    chosen = {d.name for _, d, _ in select_samples(p.samples(), only, p.settings.classification.pool)}
    return {f: v for f, v in folders.items() if v[1] in chosen}


def same_folder(a: str | Path, b: str | Path) -> bool:
    try:
        return bool(str(a).strip() and str(b).strip()) and Path(a).resolve() == Path(b).resolve()
    except OSError:
        return False


class Job:
    """One run of `imagestream-cli run-project` at a time; its output is kept for the page, and which result folders
    it writes (`folders`) since when (`started`, wall clock), so that the overview can tell old results from new."""

    def __init__(self) -> None:
        self.proc: subprocess.Popen[str] | None = None
        self.lines: list[str] = []
        self.progress = ""
        self.folder = ""
        self.folders: dict[str, tuple[str, str]] = {}
        self.started = 0.0
        self.ended: float | None = None
        self.stop_requested = False
        self.lock = threading.Lock()

    @property
    def running(self) -> bool:
        return self.proc is not None and self.proc.poll() is None

    def start(self, project_file: Path, only: list[str] | None = None) -> None:
        if self.running:
            raise RuntimeError("A run is still going")
        exe = Path(sys.executable)
        python = exe.with_name("python.exe") if exe.name.lower() == "pythonw.exe" else exe
        cmd = [str(python), "-m", "imagestream_analysis.cli", "run-project", str(project_file)]
        for name in only or []:  # only these datasets (with pools per experiment: their experiments)
            cmd += ["--dataset", name]
        try:
            folders = run_folders(project_file, only)
        except Exception as e:  # noqa: BLE001 - only for the page; the run itself reports the problem
            log.warning("Could not tell which datasets the run writes: %s", e)
            folders = {}
        with self.lock:
            self.lines, self.progress, self.folder = [], "", str(project_file.parent)
            self.folders, self.started, self.ended, self.stop_requested = folders, time.time(), None, False
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
        self.proc.wait()
        with self.lock:
            self.ended = time.time()

    def stop(self) -> None:
        if self.running and self.proc is not None:
            self.stop_requested = True
            self.proc.terminate()

    def status(self, since: int) -> dict[str, Any]:
        with self.lock:
            running = self.running
            end = time.time() if running or self.ended is None else self.ended
            return {"running": running, "lines": self.lines[since:], "n": len(self.lines),
                    "progress": self.progress, "folder": self.folder, "datasets": len(self.folders),
                    "seconds": round(end - self.started, 1) if self.started else None,
                    "stopped": self.stop_requested,
                    "returncode": None if self.proc is None else self.proc.poll()}

    def result_status(self, folder: Path, table: pd.DataFrame) -> list[str]:
        """Per row of the overview of `folder`: "running" (in the run going on, not done yet), "done" (written by
        the run going on or the last one), "old" (in the last run, which ended before writing it: stopped, or
        failed), or "" (not concerned)."""
        with self.lock:
            if not self.folders or not same_folder(folder, self.folder):
                return [""] * len(table)
            running, folders, started = self.running, self.folders, self.started
        out = []
        for rel in table["folder"]:
            rel = Path(str(rel)).as_posix()
            if rel not in folders:
                out.append("")
                continue
            try:
                fresh = (folder / rel / "summary.json").stat().st_mtime >= started
            except OSError:
                fresh = False
            out.append("done" if fresh else "running" if running else "old")
        return out

    def waiting(self, folder: Path, present: set[str]) -> list[dict[str, str]]:
        """The datasets of the run going on that have no result folder yet (a first run): rows for the overview."""
        with self.lock:
            if not self.running or not same_folder(folder, self.folder):
                return []
            return [{"experiment": e, "dataset": d, "folder": rel, "status": "waiting"}
                    for rel, (e, d) in self.folders.items() if rel not in present]


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


def apply_from_dashboard(data: dict[str, Any], job: Job) -> dict[str, Any]:
    """The Apply button of a dashboard opened from the overview: its subtype thresholds (by marker name; None = back
    to automatic) and membrane pair into the project, then the project is run again (Cellpose is skipped where its
    input and settings did not change). dry_run: only how many datasets it concerns."""
    root = SERVED.get("root")
    if root is None or not (root / PROJECT_FILE).is_file():
        raise ValueError("this dashboard does not belong to the project that is open on the import page")
    p = Project.load(root / PROJECT_FILE)
    n = len(p.samples())
    if data.get("dry_run"):
        return {"datasets": n, "running": job.running}
    if job.running:
        raise RuntimeError("a run is still going")
    thresholds = dict(p.settings.classification.thresholds)
    for name, value in (data.get("thresholds") or {}).items():
        if value is None:
            thresholds.pop(name, None)
        else:
            thresholds[name] = float(value)
    p.settings.classification.thresholds = thresholds
    if data.get("membrane"):  # the pair of one marker (the one shown in the dashboard); the others stay
        pair = data["membrane"]
        cells, partners = [str(c) for c in pair["cells"]], [str(c) for c in pair["partners"]]
        if pair.get("marker"):
            p.settings.membrane.pairs[str(pair["marker"])] = {"cells": cells, "partners": partners}
        else:  # a page from before 2026-11-03: one pair for all markers
            p.settings.membrane.cells, p.settings.membrane.partners = cells, partners
    job.start(p.save())
    return {"started": True, "datasets": n}


CELLPOSE_CACHE = "cellpose_check.json"  # in STATE_DIR: the last result per Python
_CELLPOSE_PROBE = (
    "import json\n"
    "info = {}\n"
    "try:\n"
    "    from importlib.metadata import version\n"
    "    info['cellpose'] = version('cellpose')\n"
    "    import torch\n"
    "    info['torch'] = torch.__version__\n"
    "    info['gpu'] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None\n"
    "except Exception as e:\n"
    "    info['error'] = f'{type(e).__name__}: {e}'\n"
    "print(json.dumps(info))\n"
)


def check_cellpose(setting: str | None, fresh: bool = False) -> dict[str, Any]:
    """The Cellpose version and GPU of the environment in 'Python with Cellpose' (empty: this one), as in
    FAST-HIPPOS. Starting that Python and importing PyTorch takes seconds, so the last result is returned at once
    (with "cached": true) unless fresh is set; the page then asks again with fresh=True in the background."""
    from .segmentation import cellpose_python

    try:
        exe = cellpose_python(setting)
    except ValueError as e:
        return {"error": str(e)}
    cache_file = STATE_DIR / CELLPOSE_CACHE
    try:
        cache: dict[str, Any] = json.loads(cache_file.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        cache = {}
    if not fresh and exe in cache:
        cached: dict[str, Any] = cache[exe] | {"cached": True}
        return cached
    try:
        r = subprocess.run([exe, "-I", "-c", _CELLPOSE_PROBE], capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=180, check=False, creationflags=NO_WINDOW)
        lines = r.stdout.strip().splitlines()
        result: dict[str, Any] = json.loads(lines[-1]) if lines else {"error": r.stderr[-500:]}
    except Exception as e:  # noqa: BLE001 - shown on the page
        result = {"error": f"{type(e).__name__}: {e}"}
    result["python"] = exe
    cache[exe] = result
    try:
        STATE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache, indent=1), encoding="utf-8")
    except OSError:
        pass
    return result


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
                "threshold_methods": THRESHOLD_METHODS, "pools": [str(p) for p in Pool],
                "band_centres": [str(b) for b in BandCentre]}
    if api == "prefs":
        return save_prefs({k: v for k, v in data.items() if k in ("theme",)})
    if api == "browse":
        return browse(str(data.get("kind", "zip")), str(data.get("initial", "")))
    if api == "scan":
        return scan_zip(str(data["path"]))
    if api == "cellpose":
        return check_cellpose(data.get("python") or None, bool(data.get("fresh")))
    if api == "open_project":
        p = Project.load(str(data["path"]))
        save_prefs({"project": str(p.folder)})
        return project_to_json(p)
    if api == "save":
        p = project_from_json(data["project"])
        path = p.save()
        save_prefs({"project": str(p.folder)})
        table = overview(p.folder)
        if len(table):  # the comparison follows the conditions and replicates as saved
            write_comparison(p, table)
        return {"path": str(path), "problems": p.problems()}
    if api == "run":
        p = project_from_json(data["project"])
        problems = p.problems()
        if problems:
            return {"started": False, "problems": problems}
        only = [str(n) for n in data.get("datasets") or []]
        if only:
            select_samples(p.samples(), only, p.settings.classification.pool)  # an unknown name: an error here
        path = p.save()
        save_prefs({"project": str(p.folder)})
        job.start(path, only or None)
        return {"started": True, "problems": []}
    if api == "status":
        return job.status(int(data.get("since", 0)))
    if api == "apply":
        return apply_from_dashboard(data, job)
    if api == "stop":
        job.stop()
        return {"stopped": True}
    if api == "overview":
        folder = Path(str(data["folder"]))
        if folder.is_dir():
            SERVED["root"] = folder
        found = overview(folder) if folder.is_dir() else pd.DataFrame()
        busy = job.running and same_folder(folder, job.folder)
        waiting = job.waiting(folder, {Path(str(f)).as_posix() for f in found.get("folder", [])})
        if found.empty and not waiting:
            return {"columns": [], "rows": [], "interactions": [], "comparison": None, "busy": busy}
        table = found
        # the project as on the page (saved or not): which result folders are part of it, and the comparison of
        # conditions with its labels
        proj: Project | None = project_from_json(data["project"]) if data.get("project") else (
            Project.load(folder) if (folder / PROJECT_FILE).exists() else None)
        if proj is not None:
            proj.folder = folder
        status = pd.Series(job.result_status(folder, table) if len(table) else [], index=table.index, dtype=object)
        if proj is not None and len(table):  # results of an earlier setup in this folder (removed, unticked, renamed)
            status = status.where(in_project(table, proj), "not in project")
        shown = table.assign(status=status)
        if waiting:
            shown = pd.concat([shown, pd.DataFrame(waiting)], ignore_index=True)
        table = table[status != "not in project"] if len(table) else table
        contacts = folder / INTERACTIONS_FILE
        summed = pd.read_csv(contacts, sep="	") if contacts.exists() else pd.DataFrame()
        if len(summed):
            summed = summed[summed["dataset"] == "(all)"]
            if proj is not None:  # interactions.tsv of an older version may still hold datasets left out since
                summed = summed[summed["experiment"].isin({e for e, _ in result_folders(proj).values()})]
        comparison: dict[str, Any] | None = None
        if proj is not None and len(table):
            try:
                result, notes, conditions, reference = compare(proj, table)
                comparison = {"rows": json.loads(result.to_json(orient="records", force_ascii=False)),
                              "notes": notes, "conditions": conditions, "reference": reference}
            except Exception as e:  # noqa: BLE001 - the overview itself must still show
                log.exception("Comparison failed")
                comparison = {"error": f"The comparison failed: {e}"}
        return {"columns": [c for c in shown.columns if c != "status"],
                "rows": json.loads(shown.to_json(orient="records", force_ascii=False)),
                "interactions": json.loads(summed.to_json(orient="records", force_ascii=False)) if len(summed) else [],
                "comparison": comparison, "busy": busy}
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
            try:
                self.send_response(code)
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                self.wfile.write(body)
            except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
                pass  # the page went away (reloaded or closed) before the answer

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


def main_gui() -> int:
    """The entry point of installs made before launcher.py (until the next uv sync rewrites imagestream.exe)."""
    from .launcher import main_gui as start

    return start()
