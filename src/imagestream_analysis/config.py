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
    # an IDEAS export (.zip, one folder per dataset) or an ImageJ hyperstack from ImageStreamCombiner (.tif, TCYX,
    # one frame per event)
    image: Path = Path()
    output: Path = Path("output")  # results go to <output>/<dataset or image stem>/
    max_events: int | None = None  # analyze only the first N events (macro: limitNrCells); None = all
    # IDEAS zip only: the dataset (folder in the zip; "" = the only one), the channels to load (IDEAS channel
    # numbers; the ones with a role are always loaded), their names (same order as `channels`) and the pixel size
    dataset: str = ""
    channels: list[int] = field(default_factory=list)
    channel_names: list[str] = field(default_factory=list)
    pixel_size: float | None = None  # µm; required for an IDEAS zip (the files do not hold it)


@dataclass
class ChannelSettings:
    # channel numbers (the IDEAS channel Ch<N> for a zip, the position in the stack for a Combiner .tif) or names
    brightfield: Channel = 1
    marker_of_interest: Channel = 2
    tumor: Channel = 3
    # markers of the immune cells (ClassifyCellTypes 2.8: phenotypic markers); each is normalized, and their sum is
    # the immune-cell marker of the segmentation and the tumor / immune classification (macro v1.8: one T-cell marker)
    immune: list[Channel] = field(default_factory=lambda: [4])
    # markers of immune subtypes (2.8: functional markers), each with a positive / negative threshold
    functional: list[Channel] = field(default_factory=list)
    marker_location: MarkerLocation = MarkerLocation.NUCLEAR  # of the tumor and immune-cell markers
    cell_of_interest: CellOfInterest = CellOfInterest.IMMUNE_CELL  # the cell type that expresses the marker of interest

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
        return {"brightfield": [self.brightfield], "marker_of_interest": [self.marker_of_interest],
                "tumor": [self.tumor], "immune": list(self.immune), "functional": list(self.functional)}

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
        out.marker_of_interest = number(self.marker_of_interest, "marker_of_interest")
        out.tumor = number(self.tumor, "tumor")
        out.immune = [number(c, "immune marker") for c in self.immune]
        out.functional = [number(c, "functional marker") for c in self.functional]
        return out


@dataclass
class SegmentationSettings:
    method: str = "cellpose"  # "cellpose" or "labelmap" (load labelmap_file)
    labelmap_file: Path | None = None  # a labelmap_cells.zip/.tif of an earlier run (montage or stack)
    model: str = "cyto3"  # Cellpose 3: cyto3, cyto2, nuclei, ...; Cellpose 4: cpsam; or a path to a model
    diameter: float = 30.0  # pixels; 0 = let Cellpose estimate it (Cellpose 3 only)
    flow_threshold: float = 0.5
    cellprob_threshold: float = 0.0
    use_gpu: bool = True
    # Cellpose segments one montage of all events (so that it looks like a microscopy image): a grid with this many
    # events per row. 0 = about square (default, SPEC I2); the number of events = one row, as the Fiji macro.
    montage_columns: int = 0
    # Python of an environment with Cellpose (3 or 4); None = this environment
    cellpose_python: Path | None = None
    erosion: int = 1  # shrink the cells from the background with this many pixels
    min_cell_size: int = 100  # pixels (after erosion)
    max_cell_size: int = 10000
    # cells cut by the edge of their event image are always flagged (column Edge); true = leave them out
    exclude_edge_cells: bool = False


class Pool(StrEnum):
    EXPERIMENT = "experiment"  # the datasets of one experiment (project) are classified together
    DATASET = "dataset"


@dataclass
class ClassificationSettings:
    # Thresholds of the functional markers: on log(intensity - min + 1) of the background-subtracted cell means (as
    # ClassifyCellTypes 2.8), with an ImageJ auto-threshold method: "MaxEntropy" (2.8) or "Otsu"
    threshold_method: str = "MaxEntropy"
    # a threshold set by hand: marker name -> background-subtracted intensity; overrides the automatic one
    thresholds: dict[str, float] = field(default_factory=dict)
    # a detection positive in at least this many functional markers (thresholded on all detections) is not a cell
    # but an artefact; it is left out of the tumor / immune clustering. 0 = no artefacts
    artefact_min_positive: int = 0
    # in a project: classify the datasets of each experiment together (one set of thresholds, comparable
    # subtypes), or each dataset on its own
    pool: Pool = Pool.EXPERIMENT


@dataclass
class MembraneSettings:
    # straightened membrane profiles of the cells of interest and the enrichment at the interface (optional: not
    # every ImageStream analysis needs it); the per-event interface table is always made
    analysis: bool = False
    thickness: int = 3  # pixels: width of the membrane band, and reach of the interface
    straightening_width: int = 5  # pixels: line width for straightening the membranes
    max_histogram_value: float = 20.0  # standard deviations above background; the minimum is -5
    histogram_bin_width: float = 0.25
    # cells without interface: this many random control windows each (as long as the median interface)
    control_repeats: int = 10


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
