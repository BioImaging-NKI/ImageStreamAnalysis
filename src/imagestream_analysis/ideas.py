"""Reading IDEAS exports directly (replaces ImageStreamCombiner).

An export is a zip with one folder per dataset, holding one file per object and channel: `<object>_Ch<N>.ome.tif`,
or `<prefix>_<object>_Ch<N>.ome.tif` when a file-name prefix was set on export (often the name of the data file).
The files hold no pixel size and no channel names: those come from the user. Their OME metadata (operator, source
path) is not read.

Differences from ImageStreamCombiner: the events are ordered by object number (the Combiner: zip order) and the
object numbers are kept; only the channels that are needed are loaded, and every file is read once.
"""

from __future__ import annotations

import base64
import io
import logging
import math
import re
import struct
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from .io import EventStack

log = logging.getLogger(__name__)
SUFFIX = ".ome.tif"
# the object number is the last number before _Ch; anything before it is the export's prefix (2026-11-05: an export
# with a prefix was read as empty)
FILE_RE = re.compile(r"^(?:(?P<prefix>.*)_)?(?P<object>\d+)_Ch(?P<channel>\d+)" + re.escape(SUFFIX) + "$",
                     re.IGNORECASE)
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
    """The datasets (first-level folders) of an IDEAS zip. Files in the root or deeper folders are ignored. A folder
    whose files carry different prefixes (several data files exported into one folder: their object numbers may
    clash) gives one dataset per prefix, named "<folder> - <prefix>"."""
    found: dict[str, dict[str, dict[int, dict[int, str]]]] = {}  # folder -> prefix -> object -> channel -> name
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            parts = name.split("/")
            if len(parts) != 2 or not (m := FILE_RE.match(parts[1])):
                continue
            objects = found.setdefault(parts[0], {}).setdefault(m["prefix"] or "", {})
            objects.setdefault(int(m["object"]), {})[int(m["channel"])] = name
    datasets: dict[str, DatasetInfo] = {}
    for folder, prefixes in found.items():
        for prefix, objects in prefixes.items():
            label = folder if len(prefixes) == 1 else f"{folder} - {prefix or '(no prefix)'}"
            datasets[label] = DatasetInfo(label, objects)
    return dict(sorted(datasets.items()))


_TIFF_TYPES = {3: ("H", 2), 4: ("I", 4)}  # SHORT, LONG: the value types of the tags read below
_SAMPLE_FORMATS = {(1, 8): "u1", (1, 16): "u2", (1, 32): "u4", (2, 8): "i1", (2, 16): "i2", (2, 32): "i4",
                   (3, 32): "f4", (3, 64): "f8"}


def _plain_tiff(buf: bytes) -> np.ndarray | None:
    """The first image of a plain TIFF (uncompressed, one sample per pixel, strips), read straight from its first
    IFD: what IDEAS writes, about ten times faster than tifffile (2026-11-04: reading the files took a large part
    of a run). None for anything else (tifffile then reads it)."""
    order = {b"II": "<", b"MM": ">"}.get(buf[:2])
    if order is None or len(buf) < 8 or struct.unpack_from(order + "H", buf, 2)[0] != 42:
        return None
    at = struct.unpack_from(order + "I", buf, 4)[0]
    n = struct.unpack_from(order + "H", buf, at)[0]
    tags: dict[int, tuple[int, ...]] = {}
    for tag, kind, count, value in struct.iter_unpack(order + "HHI4s", buf[at + 2 : at + 2 + 12 * n]):
        if kind not in _TIFF_TYPES:
            continue
        code, size = _TIFF_TYPES[kind]
        if size * count <= 4:
            tags[tag] = struct.unpack_from(f"{order}{count}{code}", value)
        else:
            tags[tag] = struct.unpack_from(f"{order}{count}{code}", buf, struct.unpack(order + "I", value)[0])
    try:
        width, height = tags[256][0], tags[257][0]
        if tags.get(259, (1,))[0] != 1 or tags.get(277, (1,))[0] != 1 or tags.get(284, (1,))[0] != 1:
            return None
        dtype = np.dtype(order + _SAMPLE_FORMATS[(tags.get(339, (1,))[0], tags.get(258, (1,))[0])])
        offsets, counts = tags[273], tags[279]
    except (KeyError, IndexError):
        return None
    data = buf[offsets[0] : offsets[0] + counts[0]] if len(offsets) == 1 else b"".join(
        buf[o : o + c] for o, c in zip(offsets, counts, strict=True))
    if len(data) < width * height * dtype.itemsize:
        return None
    out: np.ndarray = np.frombuffer(data, dtype=dtype, count=width * height).reshape(height, width)
    return out.astype(dtype.newbyteorder("="))


def _read(z: zipfile.ZipFile, name: str) -> np.ndarray:
    buf = z.read(name)
    fast = _plain_tiff(buf)
    if fast is not None:
        return fast
    with tifffile.TiffFile(io.BytesIO(buf)) as tif:
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
    extents = np.empty((len(objects), 4), dtype=np.int64)
    for t, cell in enumerate(images):
        h, w = cell[0].shape
        y, x = (shape[0] - h) // 2, (shape[1] - w) // 2  # as _pad
        extents[t] = y, x, y + h, x + w
        for i, im in enumerate(cell):
            data[t, i] = _pad(im, shape, round(medians[i]) if dtype.kind in "ui" else medians[i])
    names = names or [f"Ch{c}" for c in channels]
    return EventStack(data, list(names), medians, pixel_size, dataset, n_total, list(channels),
                      np.array(objects, dtype=np.int64), extents)


def combine(stacks: list[EventStack], name: str, source_names: list[str]) -> EventStack:
    """Several parts of one sample (datasets with the same channels, in the same order) as one EventStack.

    Each part keeps its own background: its pixels are shifted by (common - its own) background per channel, the
    common background being the highest of the parts. Everything downstream subtracts the background, so the values
    above background stay those of the part, and a part with a higher background is not taken as brighter. The
    images are padded again to the largest size; the events keep their order (part after part), and `sources` says
    which part each event comes from (object numbers may repeat between parts)."""
    first = stacks[0]
    for s in stacks[1:]:
        if s.channel_names != first.channel_names:
            raise ValueError(f"The parts of {name} have different channels: {first.channel_names} and "
                             f"{s.channel_names}")
        if not math.isclose(s.pixel_size, first.pixel_size):
            raise ValueError(f"The parts of {name} have different pixel sizes: {first.pixel_size} and {s.pixel_size}")
    medians = [max((s.medians or [0.0] * len(first.channel_names))[i] for s in stacks)
               for i in range(len(first.channel_names))]
    shape = (max(s.data.shape[2] for s in stacks), max(s.data.shape[3] for s in stacks))
    dtype = first.data.dtype
    integer = dtype.kind in "ui"
    top = np.iinfo(dtype).max if integer else np.inf
    total = sum(s.n_events for s in stacks)
    data = np.empty((total, len(medians), *shape), dtype=dtype)
    extents = np.empty((total, 4), dtype=np.int64)
    at = 0
    for s in stacks:
        own = s.medians or [0.0] * len(medians)
        shift = np.array([medians[i] - own[i] for i in range(len(medians))])
        if integer:
            shift = np.round(shift)
        h, w = s.data.shape[2:]
        y, x = (shape[0] - h) // 2, (shape[1] - w) // 2
        fill = [round(m) if integer else m for m in medians]
        for i in range(len(medians)):
            data[at : at + s.n_events, i] = fill[i]
            shifted = s.data[:, i].astype(np.float64 if not integer else np.int64) + shift[i]
            if integer and shifted.max(initial=0) > top:
                log.warning("%s: part %s saturates in channel %s after the background shift", name, s.name,
                            first.channel_names[i])
            data[at : at + s.n_events, i, y : y + h, x : x + w] = np.clip(shifted, 0, top).astype(dtype)
        extents[at : at + s.n_events] = s.image_extents() + np.array([y, x, y, x])
        at += s.n_events
        log.info("%s: part %s, %d events, background %s", name, s.name, s.n_events,
                 ", ".join(f"{n} {m:.4g}" for n, m in zip(first.channel_names, own, strict=True)))
    ids = [s.object_ids if s.object_ids is not None else np.arange(1, s.n_events + 1) for s in stacks]
    sources = np.concatenate([np.full(s.n_events, k, dtype=np.int64) for k, s in enumerate(stacks)])
    return EventStack(data, list(first.channel_names), medians, first.pixel_size, name,
                      sum(s.n_events_in_file for s in stacks), first.channel_numbers,
                      np.concatenate(ids).astype(np.int64), extents, sources, list(source_names))


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
