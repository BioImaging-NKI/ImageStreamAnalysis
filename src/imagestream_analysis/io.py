"""Reading the input stack (from ImageStreamCombiner) and labelmaps; writing tables and images."""

from __future__ import annotations

import io
import logging
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile

log = logging.getLogger(__name__)
RESULT_DECIMALS = 3  # as the macro (RESULT_DECIMALS)


@dataclass
class EventStack:
    """All events of one sample: data[t, c, y, x] (t = event, c = channel).

    Channels are addressed by number: the position in an ImageStreamCombiner stack (1..C), or the IDEAS channel
    number (Ch<N>) when read from an IDEAS zip, in which case only the requested channels are loaded."""

    data: np.ndarray  # (T, C, Y, X), as stored (usually uint16)
    channel_names: list[str]
    medians: list[float] | None  # background per channel
    pixel_size: float  # µm
    name: str
    n_events_in_file: int
    channel_numbers: list[int] | None = None  # number of each channel in data; None = 1..C
    object_ids: np.ndarray | None = None  # IDEAS object number of each event; None = unknown (Combiner stacks)

    @property
    def n_events(self) -> int:
        return int(self.data.shape[0])

    def index(self, c: int) -> int:
        numbers = self.channel_numbers or list(range(1, self.data.shape[1] + 1))
        if c not in numbers:
            raise ValueError(f"Channel {c} is not in {self.name} (channels {numbers})")
        return numbers.index(c)

    def channel(self, c: int) -> np.ndarray:
        """Channel c of all events, as float32 (the macro converts the stack to 32-bit)."""
        return self.data[:, self.index(c)].astype(np.float32)

    def median(self, c: int) -> float:
        return self.medians[self.index(c)] if self.medians else 0.0


def parse_medians(properties: dict[str, Any] | None, n_channels: int) -> list[float] | None:
    """The 'Medians' property of ImageStreamCombiner: '\\n819,13,44,...' (the macro reads line 2 of the Info)."""
    text = (properties or {}).get("Medians", "")
    values = [v for v in text.strip().split(",") if v.strip()]
    if len(values) != n_channels:
        return None
    return [float(v) for v in values]


def read_stack(path: str | Path, max_events: int | None = None) -> EventStack:
    path = Path(path)
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        axes = series.axes
        meta = tif.imagej_metadata or {}
        shape = dict(zip(axes, series.shape, strict=True))
        n_t = shape.get("T", 1)
        n_c = shape.get("C", 1)
        if set(axes) - {"T", "C", "Y", "X"}:
            raise ValueError(f"{path.name}: expected a TCYX hyperstack, got axes {axes}")
        n = n_t if max_events is None else min(n_t, max_events)
        # ImageJ hyperstacks store the pages in CZT order (channels fastest): page = t * C + c
        data = tif.asarray(key=range(n * n_c)).reshape(n, n_c, shape["Y"], shape["X"])
        xres = tif.pages.first.tags.get("XResolution")
        pixel_size = 1.0
        if xres is not None and xres.value[0]:
            pixel_size = xres.value[1] / xres.value[0]
    labels = meta.get("Labels") or []
    names = [str(x) for x in labels[:n_c]] if len(labels) >= n_c else [f"C{c + 1}" for c in range(n_c)]
    medians = parse_medians(meta.get("Properties"), n_c)
    if medians is None:
        log.warning("%s: no channel medians in the metadata; the background is taken as 0", path.name)
    return EventStack(data, names, medians, pixel_size, path.stem, n_t)


def _read_tif(path: Path) -> np.ndarray:
    if path.suffix.lower() == ".zip":  # ImageJ 'save as zip': one tif inside
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.lower().endswith((".tif", ".tiff")))
            return tifffile.imread(io.BytesIO(z.read(name)))
    return tifffile.imread(path)


def read_labelmap(path: str | Path, n_events: int, height: int, width: int) -> np.ndarray:
    """A labelmap as (T, Y, X) int32. Accepts a montage as given to Cellpose (a grid of events, row by row; the
    macro's single row is a grid of one row) or a stack. Extra events beyond n_events are dropped."""
    a = _read_tif(Path(path))
    if a.ndim == 2:
        if a.shape[0] % height or a.shape[1] % width:
            raise ValueError(f"Labelmap of {a.shape[1]}x{a.shape[0]} px is not a grid of {width}x{height} px events")
        a = montage_to_stack(a, height, width)
    if a.shape[0] < n_events or a.shape[1:] != (height, width):
        raise ValueError(f"Labelmap has shape {a.shape}, expected at least ({n_events}, {height}, {width})")
    return a[:n_events].astype(np.int32)


def grid_columns(n: int, height: int, width: int) -> int:
    """Columns for a montage of n events of height x width px that is about square."""
    return max(1, min(n, round(math.sqrt(n * height / width))))


def stack_to_montage(stack: np.ndarray, columns: int | None = None, fill: float = 0) -> np.ndarray:
    """(T, Y, X[, C]) -> a grid of `columns` events per row, row by row (Fiji: Make Montage); empty tiles at the end
    get `fill`. columns=None: one row (as the macro)."""
    t, h, w = stack.shape[:3]
    cols = t if columns is None else columns
    rows = -(-t // cols)
    tiles = np.full((rows * cols, *stack.shape[1:]), fill, dtype=stack.dtype)
    tiles[:t] = stack
    grid = tiles.reshape(rows, cols, h, w, *stack.shape[3:]).swapaxes(1, 2)
    return grid.reshape(rows * h, cols * w, *stack.shape[3:])


def montage_to_stack(montage: np.ndarray, height: int, width: int) -> np.ndarray:
    """A grid montage (Y, X[, C]) -> all its tiles, row by row, (rows*cols, height, width[, C])."""
    rows, cols = montage.shape[0] // height, montage.shape[1] // width
    grid = montage.reshape(rows, height, cols, width, *montage.shape[2:]).swapaxes(1, 2)
    return grid.reshape(rows * cols, height, width, *montage.shape[2:])


def save_table(df: pd.DataFrame, path: Path, index: bool = False) -> None:
    """Tab-separated, NaN written as 'NaN' and floats with 3 decimals, as the macro's tables."""
    df.to_csv(path, sep="\t", index=index, na_rep="NaN", float_format=f"%.{RESULT_DECIMALS}f", encoding="utf-8")


def save_stack(stack: np.ndarray, path: Path, pixel_size: float) -> None:
    """A (T, Y, X) stack as ImageJ tif; labels above 65535 are stored as 32-bit float (ImageJ has no int32)."""
    dtype = np.uint16 if stack.dtype.kind in "iub" and stack.max(initial=0) < 65536 else np.float32
    tifffile.imwrite(path, stack.astype(dtype), imagej=True, resolution=(1 / pixel_size, 1 / pixel_size),
                     metadata={"axes": "TYX", "unit": "um"}, compression="zlib")
