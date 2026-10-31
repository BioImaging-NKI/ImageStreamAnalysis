"""Reading IDEAS exports directly (replaces ImageStreamCombiner).

An export is a zip with one folder per dataset, holding one file per object and channel: `<object>_Ch<N>.ome.tif`.
The files hold no pixel size and no channel names: those come from the user. Their OME metadata (operator, source
path) is not read.

Differences from ImageStreamCombiner: the events are ordered by object number (the Combiner: zip order) and the
object numbers are kept; only the channels that are needed are loaded, and every file is read once.
"""

from __future__ import annotations

import base64
import io
import logging
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from .io import EventStack

log = logging.getLogger(__name__)
SUFFIX = ".ome.tif"
FILE_RE = re.compile(r"^(?P<object>\d+)_Ch(?P<channel>\d+)" + re.escape(SUFFIX) + "$", re.IGNORECASE)
logging.getLogger("tifffile").setLevel(logging.ERROR)  # IDEAS files have tags that tifffile warns about


@dataclass
class DatasetInfo:
    name: str
    files: dict[int, dict[int, str]] = field(default_factory=dict)  # object -> channel -> name in the zip

    @property
    def channels(self) -> list[int]:
        return sorted({c for chans in self.files.values() for c in chans})

    def objects(self, channels: list[int] | None = None) -> list[int]:
        """Objects that have all `channels` (default: all channels of the dataset), in numerical order."""
        need = set(channels or self.channels)
        return sorted(o for o, chans in self.files.items() if need <= chans.keys())


def scan(path: str | Path) -> dict[str, DatasetInfo]:
    """The datasets (first-level folders) of an IDEAS zip. Files in the root or deeper folders are ignored."""
    datasets: dict[str, DatasetInfo] = {}
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            parts = name.split("/")
            if len(parts) != 2 or not (m := FILE_RE.match(parts[1])):
                continue
            ds = datasets.setdefault(parts[0], DatasetInfo(parts[0]))
            ds.files.setdefault(int(m["object"]), {})[int(m["channel"])] = name
    return dict(sorted(datasets.items()))


def _read(z: zipfile.ZipFile, name: str) -> np.ndarray:
    with tifffile.TiffFile(io.BytesIO(z.read(name))) as tif:
        return tif.pages.first.asarray()


def _pad(image: np.ndarray, shape: tuple[int, int], fill: float) -> np.ndarray:
    """Centre image in an array of `shape`, filled with `fill` (as ImageStreamCombiner)."""
    out = np.full(shape, fill, dtype=image.dtype)
    h, w = image.shape
    y, x = (shape[0] - h) // 2, (shape[1] - w) // 2
    out[y : y + h, x : x + w] = image
    return out


def read_dataset(path: str | Path, dataset: str, channels: list[int], names: list[str] | None = None,
                 pixel_size: float = 1.0, max_events: int | None = None) -> EventStack:
    """The events of one dataset as an EventStack with the IDEAS channel numbers `channels`.

    Objects that miss one of these channels are left out (logged). The background of a channel is the median over
    the objects of their median, and pads the images to the size of the largest one."""
    path = Path(path)
    info = scan(path).get(dataset)
    if info is None:
        raise ValueError(f"{path.name} has no dataset '{dataset}'")
    missing = [c for c in channels if c not in info.channels]
    if missing:
        raise ValueError(f"Dataset '{dataset}' has no channel(s) {missing} (it has {info.channels})")
    objects = info.objects(channels)
    if len(objects) < len(info.files):
        log.warning("%s: %d of %d objects miss one of the channels %s and are left out", dataset,
                    len(info.files) - len(objects), len(info.files), channels)
    n_total = len(objects)
    if max_events is not None:
        objects = objects[:max_events]
    with zipfile.ZipFile(path) as z:
        images = [[_read(z, info.files[o][c]) for c in channels] for o in objects]
    if not images:
        raise ValueError(f"Dataset '{dataset}' has no objects with channels {channels}")
    medians = [float(np.median([np.median(cell[i]) for cell in images])) for i in range(len(channels))]
    shape = (max(im[0].shape[0] for im in images), max(im[0].shape[1] for im in images))
    dtype = images[0][0].dtype
    data = np.empty((len(objects), len(channels), *shape), dtype=dtype)
    for t, cell in enumerate(images):
        for i, im in enumerate(cell):
            data[t, i] = _pad(im, shape, round(medians[i]) if dtype.kind in "ui" else medians[i])
    names = names or [f"Ch{c}" for c in channels]
    return EventStack(data, list(names), medians, pixel_size, dataset, n_total, list(channels),
                      np.array(objects, dtype=np.int64))


def _png_url(image: np.ndarray) -> str:
    buf = io.BytesIO()
    Image.fromarray(image).save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def thumbnails(path: str | Path, dataset: str, n: int = 8, size: int = 64) -> dict[int, str]:
    """For each channel, a strip of `n` example objects (the first objects with all channels, autoscaled per
    channel) as a PNG data URL, to recognise the channels when naming them."""
    info = scan(path)[dataset]
    objects = info.objects()[:n] or sorted(info.files)[:n]
    out: dict[int, str] = {}
    with zipfile.ZipFile(path) as z:
        for c in info.channels:
            tiles = []
            for o in objects:
                if c not in info.files[o]:
                    continue
                im = _read(z, info.files[o][c]).astype(np.float32)
                y, x = max((im.shape[0] - size) // 2, 0), max((im.shape[1] - size) // 2, 0)
                tiles.append(_pad(im[y : y + size, x : x + size], (size, size), float(np.median(im))))
            if not tiles:
                continue
            strip = np.concatenate([np.pad(t, ((0, 0), (0, 2)), constant_values=np.nan) for t in tiles], axis=1)
            lo, hi = np.nanpercentile(strip, [0.5, 99.5])
            scaled = np.clip((strip - lo) / ((hi - lo) or 1.0) * 255, 0, 255)
            out[c] = _png_url(np.nan_to_num(scaled, nan=255).astype(np.uint8))
    return out
