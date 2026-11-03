"""Validation of the excess over both membranes (membrane.partner_excess; manual: analysis/membrane.md) by simulation.

Events with a measured cell (a disk) touching a partner cell. The marker lies in a thin membrane at the edge of one or
both cells, drawn on a 4x finer grid so that two touching membranes lie side by side as in a real contact; the image
is blurred (Gaussian, sigma 1.2 px: the two membranes are closer than the resolution) with noise added, and analysed
with the tool's own membrane functions (straightening width 9 px, reach 3 px, as at 40x). Three cases:

  a. the marker on both cells, uniform, no enrichment: the score (SD above own membrane) is positive, the excess ~0;
  b. the marker only on the measured cell, enriched 3x at the contact: both measures find it, about equally;
  c. the marker on both cells, enriched 3x on the measured side: the excess finds it.

    python tools/stats/simulate_excess.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from imagestream_analysis.config import MembraneSettings  # noqa: E402
from imagestream_analysis.membrane import EXCESS, analyze_membranes, membrane_profiles  # noqa: E402

H = W = 64
SUPER = 4  # the fine grid on which the membranes are drawn
BLUR = 1.2  # the blur in image pixels (sigma; about 0.6 um at 40x)
CASES = {"a": (1.0, 1.0, 1.0), "b": (1.0, 0.0, 3.0), "c": (1.0, 1.0, 3.0)}  # own, partner, enrichment at contact


def _disk(cy: float, cx: float, r: float, scale: int = 1) -> np.ndarray:
    yy, xx = np.mgrid[: H * scale, : W * scale]
    out: np.ndarray = (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r
    return out


def events(case: str, n: int = 30, seed: int = 0) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """(labels, image, cells) of n events: label 2t+1 the measured cell, 2t+2 its partner."""
    own, partner, enrich = CASES[case]
    rng = np.random.default_rng(seed)
    labels = np.zeros((n, H, W), np.int32)
    image = np.zeros((n, H, W))
    s = SUPER
    for t in range(n):
        r1, r2 = 13 + rng.uniform(-1.5, 1.5), 10 + rng.uniform(-1.5, 1.5)
        # drawn on a grid SUPER times finer: each membrane a line one fine pixel wide at the cell's edge, so that two
        # touching membranes lie side by side (in one pixel of the image), as in a real contact
        a = _disk(32 * s, 22 * s, r1 * s, s)
        b = _disk(32 * s, (22 + r1 + r2 - 4) * s, r2 * s, s) & ~a  # overlapping disks: a flat contact
        ring_a, ring_b = a & ~ndi.binary_erosion(a), b & ~ndi.binary_erosion(b)
        contact = ring_a & ndi.binary_dilation(b, iterations=s)
        fine = own * ring_a + partner * ring_b + own * (enrich - 1) * contact
        img = 100 * s * fine.reshape(H, s, W, s).mean(axis=(1, 3))  # a line of one fine pixel: 100 per pixel length
        labels[t][a[s // 2 :: s, s // 2 :: s]] = 2 * t + 1
        labels[t][b[s // 2 :: s, s // 2 :: s]] = 2 * t + 2
        image[t] = ndi.gaussian_filter(img, BLUR) + rng.normal(0, 3, (H, W))
    cells = pd.DataFrame({"frame": np.repeat(np.arange(n), 2), "label": np.arange(1, 2 * n + 1),
                          "Cell type": ["Tumor", "Immune"] * n})
    return labels, image, cells


def simulate(case: str, n: int = 30, seed: int = 0, projection: str = "max") -> tuple[float, float]:
    """The medians of the score (SD above own membrane) and of the excess over both membranes; projection across
    the band: "max" (the macro's) or "sum"."""
    labels, image, cells = events(case, n, seed)
    profiles = membrane_profiles(labels, cells, image, 0.0, 9, 3)
    if projection == "sum":
        profiles = profiles.summed()
    table = analyze_membranes(cells, profiles, "Tumor", MembraneSettings(control_repeats=2),
                              np.random.default_rng(seed)).table
    return float(table["stdDevs above mean"].median()), float(table[EXCESS].median())


if __name__ == "__main__":
    for case, text in [("a", "both cells, uniform, no enrichment"), ("b", "measured cell only, enriched"),
                       ("c", "both cells, enriched on the measured side")]:
        for projection in ("max", "sum"):
            score, excess = simulate(case, projection=projection)
            print(f"{case}. {text} ({projection} across the band): SD above own membrane {score:.2f}, "
                  f"excess above both membranes {excess:.2f}")
