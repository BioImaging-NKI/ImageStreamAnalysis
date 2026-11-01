"""The example data of the manual: a synthetic IDEAS export, so that the manual shows no research data.

    python tools/manual/make_example_data.py [folder]      (default: data/manual, gitignored)

Writes `example_export.zip` with two datasets ("Donor 1", "Donor 2"), one folder each, `<object>_Ch<N>.ome.tif` per
object and channel, as IDEAS exports them, and `truth.csv` (the cells that were drawn, for checks). Channels:

    Ch1 BF          brightfield: cells darker, with a dark rim
    Ch2 Marker      the marker of interest: on the membrane of the tumor cells, enriched where they touch an
                    immune cell that is positive for "Subtype A"; faint on the immune cells
    Ch3 Tumor       tumor cells (a cytoplasmic dye)
    Ch4 Immune      immune cells (a cytoplasmic dye)
    Ch7 Subtype A   on 35 % of the immune cells
    Ch8 Subtype B   on 25 % of the other immune cells (never both)
    Ch9 Unused      background only (a channel the user leaves unnamed)

Events: mostly one tumor cell with one immune cell; also single cells, two cells of one type, and three or four
cells. Subtype A cells sit next to the tumor cell more often than chance in events of three or more cells. About 2 %
of the detections are artefacts: bright in Subtype A and Subtype B. Every image has its own size, as in IDEAS.
Deterministic (fixed seeds).
"""

from __future__ import annotations

import csv
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

CHANNELS = {1: "BF", 2: "Marker", 3: "Tumor", 4: "Immune", 7: "Subtype A", 8: "Subtype B", 9: "Unused"}
BACKGROUND = {1: 820, 2: 14, 3: 40, 4: 45, 7: 30, 8: 28, 9: 25}
DATASETS = {"Donor 1": (11, 320, 1.0), "Donor 2": (12, 280, 0.85)}  # seed, objects, brightness
# event compositions (tumor, immune) and their weights
COMPOSITIONS = [((1, 1), 0.46), ((2, 1), 0.12), ((1, 2), 0.14), ((1, 0), 0.06), ((0, 1), 0.06), ((2, 0), 0.04),
                ((0, 2), 0.04), ((2, 2), 0.04), ((1, 3), 0.04)]


def ellipse(shape: tuple[int, int], cy: float, cx: float, ry: float, rx: float, angle: float) -> np.ndarray:
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    c, s = np.cos(angle), np.sin(angle)
    u, v = (xx - cx) * c + (yy - cy) * s, -(xx - cx) * s + (yy - cy) * c
    out: np.ndarray = (u / rx) ** 2 + (v / ry) ** 2 <= 1
    return out


def place(rng: np.random.Generator, kinds: list[str]) -> list[dict]:
    """Cells as touching discs: each new cell against an earlier one, without overlap. Subtype A cells are put
    against a tumor cell when there is one."""
    cells: list[dict] = []
    for kind in kinds:
        r = rng.uniform(11, 14) if kind == "tumor" else rng.uniform(6, 8)
        for _ in range(200):
            if not cells:
                cy, cx = 0.0, 0.0
            else:
                tumors = [c for c in cells if c["kind"] == "tumor"]
                prefer = tumors and (kind == "A" and rng.random() < 0.85)
                anchor = tumors[rng.integers(len(tumors))] if prefer else cells[rng.integers(len(cells))]
                a = rng.uniform(0, 2 * np.pi)
                d = anchor["r"] + r - 0.5
                cy, cx = anchor["cy"] + d * np.sin(a), anchor["cx"] + d * np.cos(a)
            if all(np.hypot(cy - c["cy"], cx - c["cx"]) >= c["r"] + r - 0.6 for c in cells):
                cells.append({"kind": kind, "cy": cy, "cx": cx, "r": r})
                break
    return cells


def draw_event(rng: np.random.Generator, cells: list[dict],
               brightness: float) -> tuple[dict[int, np.ndarray], np.ndarray]:
    margin = rng.integers(7, 15)
    ys = [c["cy"] for c in cells]
    xs = [c["cx"] for c in cells]
    rs = [c["r"] for c in cells]
    y0 = min(y - r for y, r in zip(ys, rs, strict=True)) - margin
    x0 = min(x - r for x, r in zip(xs, rs, strict=True)) - margin
    h = int(max(y + r for y, r in zip(ys, rs, strict=True)) - y0 + margin)
    w = int(max(x + r for x, r in zip(xs, rs, strict=True)) - x0 + margin)
    shape = (max(h, 40), max(w, 40))
    labels = np.zeros(shape, np.int32)
    for k, c in enumerate(cells, start=1):
        m = ellipse(shape, c["cy"] - y0, c["cx"] - x0, c["r"] * rng.uniform(0.92, 1.0), c["r"],
                    rng.uniform(0, np.pi)) & (labels == 0)
        labels[m] = k
    img = {ch: np.zeros(shape) for ch in CHANNELS}
    rim = {k: (labels == k) & ~ndi.binary_erosion(labels == k, iterations=2) for k in range(1, len(cells) + 1)}
    for k, c in enumerate(cells, start=1):
        own = labels == k
        if c["kind"] == "tumor":
            img[3][own] += rng.lognormal(np.log(420), 0.3) * brightness
            img[2][rim[k]] += rng.lognormal(np.log(60), 0.35) * brightness
            # the marker of interest gathers where the cell touches a Subtype A cell
            for j, d in enumerate(cells, start=1):
                if d["kind"] == "A" and ndi.binary_dilation(labels == j, iterations=2)[own].any():
                    zone = own & ndi.binary_dilation(labels == j, iterations=3)
                    img[2][zone] += rng.lognormal(np.log(110), 0.3) * brightness
        else:
            img[4][own] += rng.lognormal(np.log(380), 0.3) * brightness
            img[2][rim[k]] += rng.lognormal(np.log(6), 0.4) * brightness
            if c["kind"] in ("A", "artefact"):
                img[7][own] += rng.lognormal(np.log(160), 0.35) * brightness
            if c["kind"] in ("B", "artefact"):
                img[8][own] += rng.lognormal(np.log(140), 0.35) * brightness
        # brightfield: darker cell, dark rim, texture
        img[1][own] -= 120 + rng.normal(0, 25, own.sum())
        img[1][rim[k]] -= 160
    out = {}
    for ch, a in img.items():
        if ch != 1:
            a = ndi.gaussian_filter(a, 1.2)  # the optics
            a = rng.poisson(np.maximum(a, 0) + BACKGROUND[ch]).astype(float)
        else:
            a = ndi.gaussian_filter(a, 0.8) + BACKGROUND[1] + rng.normal(0, 12, shape)
        out[ch] = np.clip(a, 0, 4095).astype(np.uint16)
    return out, labels


def immune_kind(rng: np.random.Generator) -> str:
    u = rng.random()
    return "artefact" if u < 0.03 else "A" if u < 0.37 else "B" if u < 0.37 + 0.63 * 0.25 else "negative"


def main(folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    with zipfile.ZipFile(folder / "example_export.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for dataset, (seed, n, brightness) in DATASETS.items():
            rng = np.random.default_rng(seed)
            weights = np.array([w for _, w in COMPOSITIONS])
            objects = np.sort(rng.choice(np.arange(1, 6 * n), n, replace=False))  # IDEAS numbers are not contiguous
            for obj in objects:
                n_tumor, n_immune = COMPOSITIONS[rng.choice(len(COMPOSITIONS), p=weights / weights.sum())][0]
                kinds = ["tumor"] * n_tumor + [immune_kind(rng) for _ in range(n_immune)]
                cells = place(rng, kinds)
                images, _ = draw_event(rng, cells, brightness)
                for ch, a in images.items():
                    buf = io.BytesIO()
                    tifffile.imwrite(buf, a, ome=True)
                    z.writestr(f"{dataset}/{obj}_Ch{ch}.ome.tif", buf.getvalue())
                rows += [{"dataset": dataset, "object": int(obj), "kind": c["kind"]} for c in cells]
    with (folder / "truth.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["dataset", "object", "kind"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"{folder / 'example_export.zip'}: {sum(DATASETS[d][1] for d in DATASETS)} objects, {len(rows)} cells")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[2] / "data" / "manual")
