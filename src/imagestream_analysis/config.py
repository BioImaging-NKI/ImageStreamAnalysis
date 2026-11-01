"""Run settings, loaded from a TOML file (replaces the Fiji script-parameter dialog).

Conventions:
- Channels are 1-based, as in Fiji. A channel can also be given by its name (`tumor = 'Tumor'`): the channel names
  of the stack, or `input.channel_names` of an IDEAS zip. A panel by name fits zips whose channels are numbered
  differently.
- Frames (events) are 0-based indices in the code; the output tables show them 1-based, as the macro does.
- Sizes are in pixels unless the name says otherwise.
"""

from __future__ import annotations

import copy
import tomllib
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, get_args, get_origin, get_type_hints


class MarkerLocation(StrEnum):
    NUCLEAR = "nuclear"
    MEMBRANE = "membrane"


class CellOfInterest(StrEnum):
    IMMUNE_CELL = "Immune cell"
    TUMOR_CELL = "Tumor cell"

    @classmethod
    def _missing_(cls, value: object) -> CellOfInterest | None:
        return cls.IMMUNE_CELL if value == "T cell" else None  # the name in settings before 2026-10-31


Channel = int | str  # a channel number (1-based), or a channel name


@dataclass
class InputSettings:
    # an IDEAS export (.zip, one folder per dataset) or an ImageJ hyperstack from ImageStreamCombiner (.tif, one
    # frame per event)
    image: Path = Path()
    output: Path = Path("output")  # the results go to a folder named after the dataset (or stack) in this folder
    max_events: int | None = None  # analyse only the first N events (by object number); empty = all
    dataset: str = ""  # IDEAS zip: the dataset (folder in the zip); empty = the only one
    # IDEAS zip: the channels to load (IDEAS numbers Ch<N>); the channels with a role are always loaded
    channels: list[int] = field(default_factory=list)
    channel_names: list[str] = field(default_factory=list)  # IDEAS zip: their names, in the same order as channels
    pixel_size: float | None = None  # µm per pixel; required for an IDEAS zip (its files do not hold it)


@dataclass
class ChannelSettings:
    # a channel is a number (the IDEAS channel Ch<N> for a zip, the position in the stack for a Combiner .tif) or a
    # name (of input.channel_names, or of the stack)

    brightfield: Channel = 1  # the brightfield channel
    # the marker whose enrichment at the interface is measured; empty ("") = none (only cell types and contacts; the
    # membrane analysis needs one)
    marker_of_interest: Channel = 2
    tumor: Channel = 3  # the tumor-cell marker
    # the immune-cell marker(s), a list; several are each normalized and summed (2.8: phenotypic markers)
    immune: list[Channel] = field(default_factory=lambda: [4])
    # immune subtype markers, a list: positive or negative per immune cell (2.8: functional markers)
    functional: list[Channel] = field(default_factory=list)
    # where the tumor and immune markers stain: "nuclear" (nucleus or cytoplasm) or "membrane"; sets the image for
    # Cellpose
    marker_location: MarkerLocation = MarkerLocation.NUCLEAR
    # the cell type that expresses the marker of interest: "Tumor cell" or "Immune cell"
    cell_of_interest: CellOfInterest = CellOfInterest.IMMUNE_CELL

    @staticmethod
    def migrate(data: dict[str, Any]) -> dict[str, Any]:
        """Settings written before 2026-10-31 have one T-cell marker: t_cell = N."""
        if "t_cell" in data:
            data = dict(data)
            if "immune" in data:
                raise ValueError("Give either channels.immune or (the old) channels.t_cell, not both")
            data["immune"] = [data.pop("t_cell")]
        return data

    def roles(self) -> dict[str, list[Channel]]:
        """Every role with its channel(s)."""
        return {"brightfield": [self.brightfield],
                "marker_of_interest": [self.marker_of_interest] if self.has_marker_of_interest else [],
                "tumor": [self.tumor], "immune": list(self.immune), "functional": list(self.functional)}

    @property
    def has_marker_of_interest(self) -> bool:
        return self.marker_of_interest not in ("", None)

    def resolved(self, numbers: list[int], names: list[str]) -> ChannelSettings:
        """A copy with every channel given by number: names are looked up in `names` (the channels `numbers`)."""
        def number(ref: Channel, role: str) -> int:
            if isinstance(ref, int):
                return ref
            what = role.replace("_", " ")
            if not ref:
                raise ValueError(f"no channel for the {what}")
            hits = [n for n, name in zip(numbers, names, strict=True) if name == ref]
            if len(hits) != 1:
                known = ", ".join(f"'{n}'" for n in names if n)
                raise ValueError(f"{what} '{ref}': {'no' if not hits else 'more than one'} channel with this name"
                                 f" (channels: {known or 'no names'})")
            return hits[0]

        out = copy.deepcopy(self)
        out.brightfield = number(self.brightfield, "brightfield")
        if self.has_marker_of_interest:
            out.marker_of_interest = number(self.marker_of_interest, "marker_of_interest")
        out.tumor = number(self.tumor, "tumor")
        out.immune = [number(c, "immune marker") for c in self.immune]
        out.functional = [number(c, "functional marker") for c in self.functional]
        return out


@dataclass
class SegmentationSettings:
    method: str = "cellpose"  # "cellpose", or "labelmap": load the labels of labelmap_file instead
    labelmap_file: Path | None = None  # labelmap_cells.tif (or the macro's .zip) of an earlier run
    model: str = "cyto3"  # Cellpose 3: cyto3, cyto2, nuclei, ...; Cellpose 4: cpsam; or the path of a model
    diameter: float = 30.0  # the typical cell diameter in pixels; 0 = Cellpose estimates it (Cellpose 3 only)
    flow_threshold: float = 0.5  # Cellpose's flow threshold (higher: more, less certain cells)
    cellprob_threshold: float = 0.0  # Cellpose's cell probability threshold (lower: larger cells)
    use_gpu: bool = True  # run Cellpose on the GPU
    # reuse the Cellpose result of the previous run of the dataset when its input image and Cellpose settings are the
    # same (labelmap_cells.tif and .json in the dataset folder); off = always run Cellpose
    reuse_labels: bool = True
    # events per row of the montage that Cellpose segments; 0 = about square; the number of events = one row, as the
    # Fiji macro
    montage_columns: int = 0
    cellpose_python: Path | None = None  # python.exe of an environment with Cellpose 3 or 4; empty = this one
    erosion: int = 1  # shrink the cells from the background by this many pixels
    min_cell_size: int = 100  # smaller cells (in pixels, after shrinking) are removed
    max_cell_size: int = 10000  # larger cells (in pixels, after shrinking) are removed
    # leave out cells cut by the edge of their event image (they are always marked: column Edge)
    exclude_edge_cells: bool = False


# the threshold methods offered on the page (any method of imagej.THRESHOLD_BINS works, plus " strict")
THRESHOLD_METHODS = ["Two populations", "MaxEntropy", "MaxEntropy strict", "Otsu"]


class Pool(StrEnum):
    EXPERIMENT = "experiment"  # the datasets of one experiment (project) are classified together
    DATASET = "dataset"


@dataclass
class ClassificationSettings:
    # the automatic threshold of the subtype markers, on log(intensity - min + 1): "Two populations" (where a
    # positive population is as likely as the negative one, or the upper noise limit of the negative cells when
    # there is no second population), "MaxEntropy" (as 2.8), "MaxEntropy strict" (at least 2 robust SDs above the
    # negative cells) or "Otsu"; in a settings file also "Yen" or "RenyiEntropy"
    threshold_method: str = "Two populations"
    # thresholds set by hand, by marker name: { "Subtype A" = 12.5 } (background-subtracted intensity); they replace
    # the automatic ones
    thresholds: dict[str, float] = field(default_factory=dict)
    # a detection positive in at least this many subtype markers is an artefact, not a cell; empty = all of them
    # (with two or more markers); 0 = no artefacts
    artefact_min_positive: int | None = None
    # in a project: "experiment" classifies the datasets of each experiment together (one set of thresholds),
    # "dataset" each dataset on its own
    pool: Pool = Pool.EXPERIMENT


STRAIGHTENING_UM = 4.5  # the default straightening width in µm (manual: analysis/membrane.md)
MEMBRANE_UM = 1.5  # the default membrane band and touching distance in µm (manual: analysis/interactions.md)


def membrane_px(thickness: int | None, pixel_size: float) -> int:
    """The membrane band and touching distance in pixels: as set, or the number nearest MEMBRANE_UM (at least 2)."""
    return int(thickness) if thickness else max(2, round(MEMBRANE_UM / pixel_size))


def straightening_px(width: int | None, pixel_size: float) -> int:
    """The straightening width in pixels: as set, or the odd number of pixels nearest STRAIGHTENING_UM (at least 5)."""
    if width:
        return int(width)
    return max(5, 2 * round((STRAIGHTENING_UM / pixel_size - 1) / 2) + 1)


@dataclass
class MembraneSettings:
    # the membrane profiles and the enrichment at the interface (the per-event interface table is always made)
    analysis: bool = False
    # the cells whose membrane is measured, a list of: a cell type ("Tumor", "Immune"), an immune subtype as in the
    # tables ("Immune A+ B+"), or a marker ("A+": immune cells positive for A); empty = channels.cell_of_interest
    cells: list[str] = field(default_factory=list)
    # the partners whose contact is the interface, as for cells; empty = the other cell type
    partners: list[str] = field(default_factory=list)
    # pixels: the width of the membrane band, and the distance at which cells touch; empty = about 1.5 µm (2 px at
    # 20x, 3 at 40x, 5 at 60x)
    thickness: int | None = None
    # pixels: the width of the band across the outline that is straightened; empty = about 4.5 µm (5 px at 20x, 9 at
    # 40x, 13 at 60x), so that the band reaches the membrane where the outline lies a little outside it
    straightening_width: int | None = None
    max_histogram_value: float = 20.0  # the end of the histograms (SD above the mean); they start at -5
    histogram_bin_width: float = 0.25  # the bin width of the histograms (SD)
    # random control windows per cell without an interface, each as long as the median interface
    control_repeats: int = 10
    # take the profile only from the part of the band inside the measured cell (a sensitivity test, not a separation:
    # touching membranes are far closer than the resolution, so every pixel there mixes both)
    inside_only: bool = False


@dataclass
class Settings:
    input: InputSettings = field(default_factory=InputSettings)
    channels: ChannelSettings = field(default_factory=ChannelSettings)
    segmentation: SegmentationSettings = field(default_factory=SegmentationSettings)
    classification: ClassificationSettings = field(default_factory=ClassificationSettings)
    membrane: MembraneSettings = field(default_factory=MembraneSettings)

    @classmethod
    def load(cls, path: str | Path) -> Settings:
        path = Path(path)
        with path.open("rb") as f:
            data = tomllib.load(f)
        settings: Settings = _from_dict(cls, data)
        settings.resolve_paths(path.parent)
        return settings

    def resolve_paths(self, base: Path) -> None:
        """Relative paths in a settings file are relative to that file."""
        for section in (self.input, self.segmentation):
            for f in fields(section):
                value = getattr(section, f.name)
                if isinstance(value, Path) and not value.is_absolute() and str(value) not in ("", "."):
                    setattr(section, f.name, base / value)

    def to_dict(self) -> dict[str, Any]:
        plain: dict[str, Any] = _plain(asdict(self))
        return plain

    def to_toml(self) -> str:
        lines = []
        for section, values in self.to_dict().items():
            lines.append(f"[{section}]")
            lines += [f"{k} = {_toml_value(v)}" for k, v in values.items() if v is not None]
            lines.append("")
        return "\n".join(lines)


def _plain(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, StrEnum):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _toml_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, list):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    if isinstance(v, dict):  # an inline table with quoted keys
        return "{ " + ", ".join(f"{_toml_value(str(k))} = {_toml_value(x)}" for k, x in v.items()) + " }"
    text = str(v)
    if "'" in text or "\n" in text:  # a literal string cannot hold these: basic string with escapes
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'
    return "'" + text + "'"  # literal string: Windows paths need no escaping


def _convert(tp: Any, value: Any) -> Any:
    origin = get_origin(tp)
    if origin is not None and type(None) in get_args(tp):  # X | None
        if value is None:
            return None
        tp = next(a for a in get_args(tp) if a is not type(None))
    if tp is Path:
        return Path(value)
    if isinstance(tp, type) and issubclass(tp, StrEnum):
        return tp(value)
    if tp is float and isinstance(value, int):
        return float(value)
    return value


def _from_dict(cls: Any, data: dict[str, Any]) -> Any:
    hints = get_type_hints(cls)
    if hasattr(cls, "migrate"):
        data = cls.migrate(data)
    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ValueError(f"Unknown setting(s) in [{cls.__name__}]: {', '.join(sorted(unknown))}")
    kwargs = {}
    for name, value in data.items():
        tp = hints[name]
        kwargs[name] = _from_dict(tp, value) if is_dataclass(tp) else _convert(tp, value)
    return cls(**kwargs)
