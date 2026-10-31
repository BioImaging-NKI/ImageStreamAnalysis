"""Re-implementations of the ImageJ / CLIJ2 operations that the macro uses, where their exact definition matters
for reproducing its numbers. See docs/SPEC.md for the source of each definition."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy import ndimage as ndi

CROSS = ndi.generate_binary_structure(2, 1)  # CLIJ2 *Sphere* in 2-D: von Neumann neighbourhood (4-connected)


def otsu_bin(histogram: np.ndarray) -> int:
    """ImageJ AutoThresholder 'Otsu' on a 256-bin histogram: the last bin of the background class."""
    h = histogram.astype(np.float64)
    total = h.sum()
    if total == 0:
        return 0
    p = h / total
    omega = np.cumsum(p)
    mu = np.cumsum(np.arange(len(p)) * p)
    mu_t = mu[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        bcv = (mu_t * omega - mu) ** 2 / (omega * (1.0 - omega))
    bcv[~np.isfinite(bcv)] = 0.0
    return int(np.argmax(bcv)) if bcv.max() > 0 else 0


def float_to_byte(image: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """ImageJ float -> byte with the display range [lo, hi]: round((v - lo) * 255 / (hi - lo)), clipped to 0-255.
    (Verified on the test data: the other variants miss Fiji's normalization factors by 0.3-1.5 %.)"""
    scale = 255.0 / (hi - lo)
    v = np.floor((image.astype(np.float64) - lo) * scale + 0.5)
    return np.clip(v, 0, 255).astype(np.uint8)


def max_entropy_bin(histogram: np.ndarray) -> int:
    """ImageJ AutoThresholder 'MaxEntropy' (Kapur, Sahoo and Wong 1985) on a histogram: the last bin of the
    background class, the split that maximizes the sum of the entropies of both classes."""
    h = histogram.astype(np.float64)
    total = h.sum()
    if total == 0:
        return 0
    p = h / total
    p1 = np.cumsum(p)
    p2 = 1.0 - p1
    eps = 2.220446049250313e-16
    nonzero = np.flatnonzero(np.abs(p1) >= eps)
    first = int(nonzero[0]) if len(nonzero) else 0
    nonzero = np.flatnonzero(np.abs(p2[first:]) >= eps)
    last = first + int(nonzero[-1]) if len(nonzero) else len(p) - 1
    # entropy terms -p*log(p/P) = -p*log(p) + p*log(P), summed with cumulative sums
    plogp = np.where(p > 0, p * np.log(np.where(p > 0, p, 1.0)), 0.0)
    c_plogp, c_p = np.cumsum(plogp), p1
    best: float = 5e-324  # ImageJ: max_ent = Double.MIN_VALUE
    threshold = -1
    for it in range(first, last + 1):
        ent_back = 0.0 if p1[it] <= 0 else float(-(c_plogp[it] / p1[it]) + c_p[it] / p1[it] * np.log(p1[it]))
        rest_plogp, rest_p = c_plogp[-1] - c_plogp[it], 1.0 - c_p[it]
        ent_obj = 0.0 if p2[it] <= 0 else float(-(rest_plogp / p2[it]) + rest_p / p2[it] * np.log(p2[it]))
        if ent_back + ent_obj > best:
            best, threshold = ent_back + ent_obj, it
    return threshold


THRESHOLD_BINS = {"Otsu": otsu_bin, "MaxEntropy": max_entropy_bin}


def auto_threshold_dark(image: np.ndarray, method: str = "Otsu") -> float:
    """Lower threshold of setAutoThreshold("<method> dark") on a 32-bit image, as getThreshold() returns it.

    For float images ImageJ resets the display range to min-max, converts to 8-bit, thresholds the 8-bit histogram
    (objects = bins above the threshold bin) and maps the 8-bit level back: lower = min + (bin+1)/255 * (max-min)."""
    finite = image[np.isfinite(image)]
    lo, hi = float(finite.min()), float(finite.max())
    if hi <= lo:
        return lo
    hist = np.bincount(float_to_byte(finite, lo, hi).ravel(), minlength=256)
    level = min(THRESHOLD_BINS[method](hist) + 1, 255)
    return lo + level / 255.0 * (hi - lo)


def otsu_threshold_dark(image: np.ndarray) -> float:
    """setAutoThreshold("Otsu dark") on a 32-bit image (the macro's normalizeImage)."""
    return auto_threshold_dark(image, "Otsu")


def median_upper(values: Iterable[float] | np.ndarray) -> float:
    """The macro's medianOfArray: element floor(n/2) of the sorted values (the upper median for even n)."""
    s = np.sort(np.asarray(values, dtype=np.float64))
    return float(s[len(s) // 2])


def circular_kernel(radius: float) -> np.ndarray:
    """The circular kernel of ImageJ RankFilters (r2 = radius^2 + 1; half-width of row dy = int(sqrt(r2 - dy^2)))."""
    r2 = int(radius * radius) + 1
    k = int(np.sqrt(r2 + 1e-10))
    kernel = np.zeros((2 * k + 1, 2 * k + 1), dtype=bool)
    for dy in range(-k, k + 1):
        dx = int(np.sqrt(r2 - dy * dy + 1e-10))
        kernel[dy + k, k - dx : k + dx + 1] = True
    return kernel


def variance_filter(image: np.ndarray, radius: float = 2) -> np.ndarray:
    """ImageJ Process > Filters > Variance (population variance over the circular kernel; edges padded with the
    nearest pixels). image: (Y, X) or a stack (T, Y, X), filtered per slice."""
    kernel = circular_kernel(radius).astype(np.float64)
    kernel /= kernel.sum()
    if image.ndim == 3:
        kernel = kernel[None]
    x = image.astype(np.float64)
    mean = ndi.correlate(x, kernel, mode="nearest")
    mean2 = ndi.correlate(x * x, kernel, mode="nearest")
    return np.maximum(mean2 - mean * mean, 0).astype(np.float32)


def erode(mask: np.ndarray, iterations: int) -> np.ndarray:
    """CLIJ2 erodeSphereSliceBySlice, repeated: 4-connected binary erosion per slice; pixels outside the image
    count as the nearest edge pixel (so the image border itself does not erode)."""
    out = mask.astype(bool)
    structure = CROSS[None] if out.ndim == 3 else CROSS
    for _ in range(iterations):
        out = ndi.binary_erosion(out, structure=structure, border_value=1)
    return out


def dilate(mask: np.ndarray, iterations: int) -> np.ndarray:
    """CLIJ2 dilateSphereSliceBySlice, repeated: 4-connected binary dilation per slice."""
    out = mask.astype(bool)
    structure = CROSS[None] if out.ndim == 3 else CROSS
    for _ in range(iterations):
        out = ndi.binary_dilation(out, structure=structure, border_value=0)
    return out


def close_index_gaps(labels: np.ndarray) -> np.ndarray:
    """CLIJ2 closeIndexGapsInLabelMap: renumber the labels 1..N, keeping their order."""
    values = np.unique(labels)
    values = values[values != 0]
    lut = np.zeros(int(labels.max(initial=0)) + 1, dtype=np.int32)
    lut[values] = np.arange(1, len(values) + 1, dtype=np.int32)
    out: np.ndarray = lut[labels]
    return out
