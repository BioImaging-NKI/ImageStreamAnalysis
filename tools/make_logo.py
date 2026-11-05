"""The logo of ImageStream Analysis: an SVG (pages, manual) and a Windows icon (shortcuts), from one drawing.

Two cells in the flow channel of the cytometer: a small immune cell (blue) touches a tumor cell (orange), each drawn
with its membrane, and their contact glows (aqua, the interface colour of the dashboard).

    python tools/make_logo.py

Writes src/imagestream_analysis/templates/logo.svg and logo.ico, and manual/images/logo.svg. The icon is the SVG
rendered by Edge (Playwright), so that both look the same.
"""

from __future__ import annotations

import io
import shutil
from pathlib import Path

from PIL import Image
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "imagestream_analysis" / "templates"

SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
<defs>
<linearGradient id="bg" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#25324a"/><stop offset="1" stop-color="#121722"/></linearGradient>
<radialGradient id="tumor" cx=".38" cy=".34" r=".75"><stop offset="0" stop-color="#ff9a68"/><stop offset=".7" stop-color="#e5602b"/><stop offset="1" stop-color="#b8441a"/></radialGradient>
<radialGradient id="immune" cx=".38" cy=".34" r=".75"><stop offset="0" stop-color="#7db6ff"/><stop offset=".7" stop-color="#3577d6"/><stop offset="1" stop-color="#20539e"/></radialGradient>
<filter id="glow" x="-50%" y="-50%" width="200%" height="200%"><feGaussianBlur stdDeviation="1.6"/></filter>
</defs>
<rect width="64" height="64" rx="14" fill="url(#bg)"/>
<path d="M7 14.5H57M7 49.5H57" stroke="#ffffff" stroke-opacity=".13" stroke-width="1.3" stroke-linecap="round"/>
<path d="M4.5 28.5H9M3 32H8.5M4.5 35.5H9" stroke="#9cc4ff" stroke-opacity=".45" stroke-width="1.4" stroke-linecap="round"/>
<circle cx="19" cy="32" r="8.6" fill="url(#immune)" stroke="#b9d7ff" stroke-width="1.6"/>
<circle cx="40.6" cy="32" r="13" fill="url(#tumor)" stroke="#ffc2a1" stroke-width="1.8"/>
<path d="M28.3 25.4A13 13 0 0 0 28.3 38.6" fill="none" stroke="#3ff2b4" stroke-width="5" stroke-linecap="round" opacity=".55" filter="url(#glow)"/>
<path d="M28.3 25.4A13 13 0 0 0 28.3 38.6" fill="none" stroke="#5cf7c4" stroke-width="2.6" stroke-linecap="round"/>
</svg>
"""  # noqa: E501 - the drawing, one element per line


def render(size: int) -> Image.Image:
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="msedge")
        page = browser.new_page(viewport={"width": size, "height": size})
        page.set_content(f'<html><body style="margin:0;background:transparent">'
                         f'{SVG.replace("<svg ", f"<svg width={size} height={size} ", 1)}</body></html>')
        shot = page.screenshot(omit_background=True)
        browser.close()
    return Image.open(io.BytesIO(shot)).convert("RGBA")


def main() -> None:
    (TEMPLATES / "logo.svg").write_text(SVG, encoding="utf-8")
    big = render(256)
    big.save(TEMPLATES / "logo.ico", sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
    shutil.copyfile(TEMPLATES / "logo.svg", ROOT / "manual" / "images" / "logo.svg")
    print(f"{TEMPLATES / 'logo.svg'}, {TEMPLATES / 'logo.ico'}, manual/images/logo.svg")


if __name__ == "__main__":
    main()
