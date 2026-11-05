"""Start the import and run page (`imagestream`, `imagestream-cli gui`) with the browser open at once.

Importing the analysis (pandas, scikit-learn, ...) takes several seconds, up to 20 s on a cold start. So this module
uses only the standard library: it opens the server and the browser first, serves a waiting page with a spinner, and
only then imports the page server (gui.py). The waiting page polls /ping and reloads itself once the real page is
there. A second start shows the running page again.
"""

from __future__ import annotations

import base64
import html
import json
import os
import secrets
import threading
import time
import traceback
import urllib.request
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

STATE_DIR = Path(os.environ.get("ISA_HOME") or Path.home() / ".imagestream_analysis")  # ISA_HOME: a test instance
INSTANCE = STATE_DIR / "gui_instance.json"  # port of a running page server, so that a second start reuses it
IDLE_TIMEOUT = 15 * 60  # s after the last request; the page pings every minute
READY = b"imagestream-analysis"  # /ping of the real page server
LOADING = b"imagestream-analysis-loading"  # /ping while it is still starting

WAITING_PAGE = b"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ImageStream Analysis</title>
<link rel="icon" href="__ISA_LOGO__">
<meta name="viewport" content="width=device-width, initial-scale=1"><style>
:root { color-scheme: light dark; }
body { margin: 0; height: 100vh; display: grid; place-items: center; font: 15px system-ui, sans-serif;
  background: Canvas; color: CanvasText; }
.box { text-align: center; }
.logo { display: block; width: 96px; height: 96px; margin: 0 auto 22px; border-radius: 22px; }
.spin { width: 44px; height: 44px; margin: 0 auto 18px; border-radius: 50%; border: 4px solid #2a78d633;
  border-top-color: #2a78d6; animation: turn .9s linear infinite; }
@keyframes turn { to { transform: rotate(360deg); } }
.sub { opacity: .65; font-size: 13px; margin-top: 6px; }
</style></head><body><div class="box"><img class="logo" src="__ISA_LOGO__" alt="">
<div class="spin"></div><div>Starting ImageStream Analysis&hellip;</div>
<div class="sub">Loading the analysis modules; this takes a few seconds.</div></div><script>
async function poll() {
  try {
    const r = await fetch("/ping", { cache: "no-store" });
    if ((await r.text()) === "imagestream-analysis") { location.reload(); return; }
  } catch (e) { /* not answering yet */ }
  setTimeout(poll, 300);
}
poll();
</script></body></html>"""
try:  # the logo (templates/logo.svg), inline; standard library only, as the rest of this module
    _LOGO = (Path(__file__).parent / "templates" / "logo.svg").read_bytes()
    WAITING_PAGE = WAITING_PAGE.replace(b"__ISA_LOGO__", b"data:image/svg+xml;base64," + base64.b64encode(_LOGO))
except OSError:
    WAITING_PAGE = WAITING_PAGE.replace(b"__ISA_LOGO__", b"data:,")


class Waiting(BaseHTTPRequestHandler):
    """The handler while the page server is being imported: the waiting page for every GET."""

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002 - quiet console
        pass

    def do_GET(self) -> None:  # noqa: N802
        body, ctype = (LOADING, "text/plain") if self.path == "/ping" else (WAITING_PAGE, "text/html; charset=utf-8")
        try:
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            pass

    def do_POST(self) -> None:  # noqa: N802
        self.send_error(503, "starting")


class Failed(Waiting):
    """The handler when the page server could not be imported: the error, instead of waiting forever."""

    message = ""

    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/ping":
            body, ctype = READY, "text/plain"  # the waiting page reloads and gets the error
        else:
            body = (b"<!doctype html><meta charset='utf-8'><title>ImageStream Analysis</title>"
                    b"<body style='font: 14px system-ui, sans-serif; margin: 24px'>"
                    b"<h2>ImageStream Analysis could not start</h2><pre style='white-space: pre-wrap'>"
                    + html.escape(self.message).encode("utf-8") + b"</pre></body>")
            ctype = "text/html; charset=utf-8"
        try:
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)
        except (ConnectionAbortedError, ConnectionResetError, BrokenPipeError):
            pass


def running_instance() -> str | None:
    try:
        url = json.loads(INSTANCE.read_text(encoding="utf-8"))["url"]
        with urllib.request.urlopen(url.rstrip("/") + "/ping", timeout=1) as r:
            return str(url) if r.read() in (READY, LOADING) else None
    except Exception:  # noqa: BLE001 - no (live) instance
        return None


def serve(port: int = 0, open_browser: bool = True) -> None:
    url = running_instance()
    if url and port == 0:  # a page server is running: show it again instead of starting a second one
        if open_browser:
            webbrowser.open(url)
        print(f"Already running: {url}")
        return
    server = ThreadingHTTPServer(("127.0.0.1", port), Waiting)
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    INSTANCE.write_text(json.dumps({"url": url}), encoding="utf-8")
    serving = threading.Thread(target=server.serve_forever, daemon=True)
    serving.start()
    try:
        if open_browser:
            webbrowser.open(url)

        try:
            from . import gui  # the slow part, while the browser shows the waiting page
        except Exception:
            Failed.message = traceback.format_exc()
            server.RequestHandlerClass = Failed  # the waiting page reloads and shows the error
            print(Failed.message)
            serving.join(600)
            server.shutdown()
            return

        job = gui.Job()
        last_seen = [time.monotonic()]
        server.RequestHandlerClass = gui.make_handler(secrets.token_urlsafe(24), job, last_seen)
        print(f"ImageStream analysis: {url} (stops {IDLE_TIMEOUT // 60} min after the page was closed, or with Quit)")

        def watchdog() -> None:  # stop when idle (the page pings every minute), never during a run
            while serving.is_alive():
                time.sleep(30)
                if not job.running and time.monotonic() - last_seen[0] > IDLE_TIMEOUT:
                    server.shutdown()
                    return

        threading.Thread(target=watchdog, daemon=True).start()
        while serving.is_alive():  # until Quit or the watchdog shuts the server down
            serving.join(1)
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
