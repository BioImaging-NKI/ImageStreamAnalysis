"""Screenshots of the user manual (manual/images), from the example project (make_example_project.py), with the
installed Microsoft Edge driven by Playwright. The import page runs in this process with a state folder of its own,
so the settings of the user are not touched.

    python tools/manual/screenshots.py [only these shots, e.g. gui_channels dash_markers]
"""

from __future__ import annotations

import json
import sys
import tempfile
import threading
import time
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import quote

from PIL import Image
from playwright.sync_api import Page, sync_playwright

from imagestream_analysis import gui

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "data" / "manual" / "project"
DATASET = "Example/Donor 1"
OUT = ROOT / "manual" / "images"
VIEW = {"width": 1440, "height": 900}
SCALE = 1.5  # device pixels per CSS pixel: sharp text, also on high-resolution screens


def save(target, name: str, want: set[str], **kw) -> None:
    if want and name not in want:
        return
    path = OUT / f"{name}.png"
    target.screenshot(path=path, animations="disabled", **kw)
    im = Image.open(path)
    im.save(path, optimize=True)
    print(f"{name}.png {im.size[0]}x{im.size[1]}")


def start_page_server(state: Path) -> tuple[str, str]:
    """The import page on a free port, with its own state folder that opens the example project."""
    gui.STATE_DIR = state
    gui.PREFS = state / "gui_prefs.json"
    state.mkdir(parents=True, exist_ok=True)
    gui.PREFS.write_text(json.dumps({"project": str(PROJECT), "theme": "light"}), encoding="utf-8")
    token = "manual"
    server = ThreadingHTTPServer(("127.0.0.1", 0), gui.make_handler(token, gui.Job(), [time.time()]))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return f"http://127.0.0.1:{server.server_address[1]}", token


def import_page_shots(page: Page, url: str, want: set[str]) -> None:
    page.goto(url + "/")
    page.wait_for_selector(".zip img.thumb")
    page.wait_for_selector("#overview table")
    page.wait_for_timeout(500)
    # a public manual shows no folder of this computer
    page.evaluate("""() => {
        for (const e of document.querySelectorAll('.zip .path')) e.textContent = 'D:\\\\data\\\\example_export.zip';
        document.getElementById('folder').value = 'D:\\\\data\\\\example project';
        for (const i of document.querySelectorAll('#settings input[type=text]'))
            if (/python/i.test(i.value)) i.value = 'D:\\\\Cellpose\\\\.venv\\\\Scripts\\\\python.exe'; }""")
    zips = page.locator("section", has=page.locator("h2", has_text="IDEAS exports"))
    save(zips, "gui_exports", want)
    channels = page.locator(".zip table").nth(1)
    save(channels, "gui_channels", want)
    save(page.locator("section", has=page.locator("h2", has_text="Analysis settings")), "gui_settings", want)
    save(page.locator("section", has=page.locator("h2", has_text="Overview")), "gui_overview", want)


def dashboard_shots(page: Page, url: str, token: str, want: set[str]) -> None:
    page.goto(f"{url}/files/{token}/{quote(DATASET)}/dashboard.html")
    page.wait_for_function("typeof TYPE !== 'undefined' && document.querySelectorAll('.tile canvas').length > 0")
    page.wait_for_timeout(1500)  # the image sheets
    page.add_style_tag(content="header, .filters { position: static !important; }")  # not over the cards in a shot
    save(page, "dash_overview", want, clip={"x": 0, "y": 0, "width": VIEW["width"], "height": 760})
    save(page.locator("#markers-card"), "dash_markers", want)
    page.click("#markers-card h2 .hq")  # the context help of a card
    page.wait_for_timeout(500)
    save(page, "dash_help", want, clip={"x": VIEW["width"] - 900, "y": 0, "width": 900, "height": 560})
    page.keyboard.press("Escape")
    page.wait_for_timeout(300)
    save(page.locator(".card", has=page.locator("#classes")), "dash_classes", want)
    save(page.locator("#contacts-card"), "dash_interactions", want)
    page.locator("#membrane-card").scroll_into_view_if_needed()
    page.wait_for_timeout(300)
    save(page.locator("#membrane-card"), "dash_membrane", want)
    # a cell with an interface, selected
    page.evaluate("""() => { const k = M.sd.findIndex(v => v != null && v > 3); select(M.cell[k]); }""")
    page.wait_for_timeout(800)
    save(page.locator("#cellcard"), "dash_cell", want)
    # the events with a contact between tumor cells and Subtype A cells
    page.evaluate("""() => { select(-1); document.getElementById('f-contacts').click();
                             setContact(['Tumor', 'Immune Subtype A+']); }""")
    page.wait_for_timeout(1500)
    gallery = page.locator("main > .card").first
    gallery.scroll_into_view_if_needed()
    page.wait_for_timeout(500)
    box = gallery.bounding_box()
    save(page, "dash_contact_filter", want, clip={**box, "height": min(box["height"], 640)})


def main(want: set[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="isa_manual_") as tmp:
        url, token = start_page_server(Path(tmp) / "state")
        with sync_playwright() as p:
            browser = p.chromium.launch(channel="msedge")
            page = browser.new_page(viewport=VIEW, device_scale_factor=SCALE, color_scheme="light")
            import_page_shots(page, url, want)
            dashboard_shots(page, url, token, want)
            browser.close()


if __name__ == "__main__":
    main(set(sys.argv[1:]))
