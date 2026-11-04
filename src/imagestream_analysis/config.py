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
import logging
import tomllib
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, get_args, get_origin, get_type_hints

log = logging.getLogger(__name__)

class MarkerLocation(StrEnum):
    NUCLEAR = "nuclear"
    MEMBRANE = "membrane"


Channel = int | str  # a channel number (1-based), or a channel name


@dataclass
class InputSettings:
    # The input file: an IDEAS export (.zip, with one folder per dataset) or an ImageJ hyperstack made by
    # ImageStreamCombiner (.tif, with one frame per event).
    image: Path = Path()
    # The results are written to a folder inside this folder, named after the dataset (or the stack).
    output: Path = Path("output")
    # Analyze only the first N events (in order of object number). Empty means all events.
    max_events: int | None = None
    # IDEAS zip only: the dataset to analyze (a folder in the zip). Empty means the only dataset in the zip.
    dataset: str = ""
    # IDEAS zip only: the channels to load, as IDEAS channel numbers (Ch<N>). Channels that have a role are always
    # loaded.
    channels: list[int] = field(default_factory=list)
    # IDEAS zip only: the names of these channels, in the same order as in channels.
    channel_names: list[str] = field(default_factory=list)
    # The pixel size in µm per pixel. It is required for an IDEAS zip, because the files in the zip do not record it.
    pixel_size: float | None = None


@dataclass
class ChannelSettings:
    # a channel is a number (the IDEAS channel Ch<N> for a zip, the position in the stack for a Combiner .tif) or a
    # name (of input.channel_names, or of the stack)

    brightfield: Channel = 1  # The brightfield channel.
    # The markers whose enrichment at the interface is measured, as a list. Every marker is analyzed in the same
    # way, and its columns and files carry its name. Empty ([]) means no marker of interest: the tool then
    # determines only the cell types and their contacts. The membrane analysis needs at least one marker of
    # interest.
    marker_of_interest: list[Channel] = field(default_factory=lambda: [2])
    tumor: Channel = 3  # The tumor-cell marker.
    # The immune-cell marker or markers, as a list. With several markers, each is normalized and the results are
    # summed (the phenotypic markers of macro 2.8).
    immune: list[Channel] = field(default_factory=lambda: [4])
    # The immune subtype markers, as a list. Each immune cell is called positive or negative for each of them (the
    # functional markers of macro 2.8).
    functional: list[Channel] = field(default_factory=list)
    # Where the tumor and immune markers stain: "nuclear" (the nucleus or the cytoplasm) or "membrane". This
    # determines the image that Cellpose segments.
    marker_location: MarkerLocation = MarkerLocation.NUCLEAR

    @staticmethod
    def migrate(data: dict[str, Any]) -> dict[str, Any]:
        """Older settings: one T-cell marker `t_cell = N` (before 2026-10-31); one marker of interest, `""` for none,
        and `cell_of_interest` (before 2026-11-03: the measured cells now follow from the data, so it is dropped)."""
        data = dict(data)
        if "t_cell" in data:
            if "immune" in data:
                raise ValueError("Give either channels.immune or (the old) channels.t_cell, not both")
            data["immune"] = [data.pop("t_cell")]
        moi = data.get("marker_of_interest")
        if moi is not None and not isinstance(moi, list):
            data["marker_of_interest"] = [] if moi == "" else [moi]
        if "cell_of_interest" in data:
            log.info("channels.cell_of_interest ('%s') is no longer used and is ignored: the measured cells are the "
                     "type that carries the marker of interest, or membrane.cells", data.pop("cell_of_interest"))
        return data

    def roles(self) -> dict[str, list[Channel]]:
        """Every role with its channel(s)."""
        return {"brightfield": [self.brightfield], "marker_of_interest": list(self.marker_of_interest),
                "tumor": [self.tumor], "immune": list(self.immune), "functional": list(self.functional)}

    @property
    def has_marker_of_interest(self) -> bool:
        return bool(self.marker_of_interest)

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
        out.marker_of_interest = list(dict.fromkeys(number(c, "marker of interest") for c in self.marker_of_interest))
        out.tumor = number(self.tumor, "tumor")
        # The immune markers are summed: a zip needs at least one of them, not all (another zip of the project may
        # name its immune marker differently).
        present = [c for c in self.immune if isinstance(c, int) or c in names]
        if self.immune and not present:
            number(self.immune[0], "immune marker")  # raises: no channel with this name
        out.immune = [number(c, "immune marker") for c in present]
        out.functional = [number(c, "functional marker") for c in self.functional]
        return out


@dataclass
class SegmentationSettings:
    # "cellpose" runs Cellpose. "labelmap" loads the labels from labelmap_file instead.
    method: str = "cellpose"
    labelmap_file: Path | None = None  # The labelmap_cells.tif (or the Fiji macro's .zip) of an earlier run.
    # The Cellpose model: cyto3, cyto2, nuclei, and so on for Cellpose 3, or cpsam for Cellpose 4. It can also be
    # the path of a model file.
    model: str = "cyto3"
    # The typical cell diameter in pixels. 0 means that Cellpose estimates it (Cellpose 3 only).
    diameter: float = 30.0
    # The flow threshold of Cellpose. A higher value gives more cells, some of them less certain.
    flow_threshold: float = 0.5
    # The cell probability threshold of Cellpose. A lower value gives larger cells.
    cellprob_threshold: float = 0.0
    use_gpu: bool = True  # Run Cellpose on the GPU.
    # Reuse the Cellpose result of the previous run of the dataset if its input image and Cellpose settings have not
    # changed (labelmap_cells.tif and its .json file in the dataset folder). When off, Cellpose always runs.
    reuse_labels: bool = True
    # The number of events per row in the montage that Cellpose segments. 0 makes the montage about square. Setting
    # it to the number of events gives a single row, as in the Fiji macro.
    montage_columns: int = 0
    # The python.exe of an environment with Cellpose 3 or 4, or the folder of that environment (a venv, conda or
    # pixi environment, or a project folder with .venv or .pixi). Empty means the environment of this program.
    cellpose_python: Path | None = None
    erosion: int = 1  # Shrink the cells by this many pixels at their border with the background.
    # Cells smaller than this (in pixels, after shrinking) are removed.
    min_cell_size: int = 100
    # Cells larger than this (in pixels, after shrinking) are removed.
    max_cell_size: int = 10000
    # Leave out cells that are cut by the edge of their event image. Such cells are always flagged in the column
    # Edge, whether or not they are left out.
    exclude_edge_cells: bool = False


# the threshold methods offered on the page (any method of imagej.THRESHOLD_BINS works, plus " strict")
THRESHOLD_METHODS = ["Two populations", "MaxEntropy", "MaxEntropy strict", "Otsu"]


class BandCentre(StrEnum):
    MEMBRANE = "membrane"  # where the membrane is found in the data, per cell type, free and at contacts
    OUTLINE = "outline"  # on the outline of the cell (macro 1.8)
    CELL = "cell"  # as membrane, but each cell's free membrane on its own peak (2026-11-04)


class Pool(StrEnum):
    EXPERIMENT = "experiment"  # the datasets of one experiment (project) are classified together
    DATASET = "dataset"


@dataclass
class ClassificationSettings:
    # The method for the automatic threshold of the subtype markers, applied to log(intensity - min + 1).
    # "Two populations" puts the threshold where a cell is as likely to belong to the positive population as to the
    # negative one; when there is no second population, it uses the upper noise limit of the negative cells.
    # "MaxEntropy" works as in macro 2.8. "MaxEntropy strict" raises the MaxEntropy threshold, if needed, to at
    # least 2 robust SDs above the negative cells. "Otsu" is the fourth choice. A settings file also accepts "Yen"
    # and "RenyiEntropy".
    threshold_method: str = "Two populations"
    # Thresholds set by hand, by marker name, for example { "Subtype A" = 12.5 }, as background-subtracted
    # intensities. They replace the automatic thresholds.
    thresholds: dict[str, float] = field(default_factory=dict)
    # A detection that is positive in at least this many subtype markers is an artefact, not a cell. Empty means
    # all subtype markers (only with two or more of them). 0 means that no detection is an artefact.
    artefact_min_positive: int | None = None
    # In a project, "experiment" classifies the datasets of each experiment together, with one set of thresholds.
    # "dataset" classifies each dataset on its own.
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
    # Measure the membrane profiles and the enrichment at the interface. The interface table per event is always
    # made, whether or not this is on.
    analysis: bool = False
    # The cells whose membrane is measured, as a list. Each item is a cell type ("Tumor", "Immune"), an immune
    # subtype as written in the tables ("Immune A+ B+"), or a marker ("A+" means the immune cells positive for A).
    # Empty means the cell type that carries the marker of interest, decided per marker from the data (the median
    # intensity of the membrane that touches no cell). Set here, it applies to every marker of interest without an
    # entry of its own in pairs (the import page writes pairs only).
    cells: list[str] = field(default_factory=list)
    # The partners whose contact forms the interface, written as for cells. Empty means the other cell type.
    partners: list[str] = field(default_factory=list)
    # The pair of each marker of interest, by marker name, for example { Marker = { cells = ['Tumor'], partners =
    # ['B+'] } }. A side that an entry leaves out comes from cells or partners above, and else is automatic.
    pairs: dict[str, dict[str, list[str]]] = field(default_factory=dict)

    @staticmethod
    def migrate(data: dict[str, Any]) -> dict[str, Any]:
        """Check membrane.pairs: { marker = { cells = [...], partners = [...] } }; a single string is a list of one."""
        pairs = data.get("pairs")
        if pairs is None:
            return data
        if not isinstance(pairs, dict):
            raise ValueError("membrane.pairs: a table by marker name, e.g. { Marker = { cells = ['Tumor'] } }")
        clean: dict[str, dict[str, list[str]]] = {}
        for marker, entry in pairs.items():
            if not isinstance(entry, dict) or set(entry) - {"cells", "partners"}:
                raise ValueError(f"membrane.pairs.{marker}: give cells and/or partners, e.g. {{ cells = ['Tumor'] }}")
            sides = {k: [v] if isinstance(v, str) else [str(x) for x in v] for k, v in entry.items()}
            clean[str(marker)] = {k: v for k, v in sides.items() if v}
        return dict(data) | {"pairs": {k: v for k, v in clean.items() if v}}

    def pair_settings(self, marker: str) -> tuple[list[str], list[str]]:
        """The cells and partners set for one marker: its entry in pairs, else cells / partners (each side on its
        own); empty: automatic (membrane.pair_of)."""
        entry = self.pairs.get(marker, {})
        return list(entry.get("cells") or self.cells), list(entry.get("partners") or self.partners)
    # In pixels: the width of the membrane band, which is also the distance at which cells count as touching. Empty
    # means about 1.5 µm (2 px at 20x, 3 px at 40x, 5 px at 60x).
    thickness: int | None = None
    # In pixels: the width of the band across the outline that is straightened. Empty means about 4.5 µm (5 px at
    # 20x, 9 px at 40x, 13 px at 60x). This width lets the band reach the membrane where the outline lies a little
    # outside it.
    straightening_width: int | None = None
    # Where the band across the membrane is centred. "membrane": where the membrane lies, found in the data: per cell
    # type, the marker of interest averaged across the outline, separately at the free membrane and at the contacts
    # with each cell type. "cell": the same, but every cell's free membrane on the peak of its own (cells differ from
    # each other), its contacts moved by as much plus the type's difference between contact and free membrane; cells
    # with too little free membrane, or without a clear peak, as "membrane". "outline": on the outline, as macro
    # 1.8. The outline lies outside a free membrane (about 1 µm) but on the contact line where cells touch, so a band
    # on the outline reaches less of the free membrane than of an interface, and an interface can look enriched
    # without cause.
    band_centre: BandCentre = BandCentre.MEMBRANE
    # The upper end of the histograms, in SD above the mean. The histograms start at -5. Empty means automatic: the
    # 99th percentile of the interface values or of the control windows, whichever is higher (at least 5), rounded up
    # to a whole number of bins. Values above the upper end are not counted.
    max_histogram_value: float | None = None
    # The bin width of the histograms, in SD. Empty means automatic: the smallest of 0.1, 0.2, 0.25, 0.5, 1, 2.5, 5,
    # 10, ... that gives at most 60 bins (usually 25 to 60).
    histogram_bin_width: float | None = None
    # The number of random control windows per cell without an interface. Each window is as long as the median
    # interface.
    control_repeats: int = 10
    # Take the profile only from the part of the band that lies inside the measured cell. This is a sensitivity
    # test, not a way to separate the two membranes: touching membranes are much closer together than the
    # resolution, so every pixel at the interface contains signal from both.
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
            # a table of tables (membrane.pairs) as sub-tables, [membrane.pairs.'Marker'], so that a pasted entry fits
            nested = {k: v for k, v in values.items() if isinstance(v, dict) and v
                      and all(isinstance(x, dict) for x in v.values())}
            lines += [f"{k} = {_toml_value(v)}" for k, v in values.items()
                      if v is not None and k not in nested and not (isinstance(v, dict) and not v and k == "pairs")]
            lines.append("")
            for k, table in nested.items():
                for name, entry in table.items():
                    lines.append(f"[{section}.{k}.{_toml_value(str(name))}]")
                    lines += [f"{kk} = {_toml_value(vv)}" for kk, vv in entry.items()]
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
