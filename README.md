# ImageStream Analysis (Python)

Analysis of tumor–immune cell interactions in imaging flow cytometry (ImageStream) data: segmentation of the cells
in every event image, classification into tumor cells and immune cells, and quantification of the enrichment of a marker
of interest at the interface between the cells.

This is the Python version of the Fiji macro
`Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm`
([BioImaging-NKI/ImageStreamAnalysis](https://github.com/BioImaging-NKI/ImageStreamAnalysis)), which is kept in
[`reference/fiji_v1.8`](reference/fiji_v1.8). Its behaviour, and every deliberate difference from it, is described in
[`docs/SPEC.md`](docs/SPEC.md).

> **Status: early development.** The analysis of macro v1.8 is ported, including the membrane profiles
> (optional), with an interactive dashboard per dataset.

## Installation

With [uv](https://docs.astral.sh/uv/):

```powershell
uv sync
```

Cellpose can run from another environment (recommended, see `segmentation.cellpose_python` in the settings), or be
installed here with `uv sync --extra cellpose3` or `uv sync --extra cellpose4`.

## Usage

The input is an IDEAS export (a zip with one folder per dataset), or an ImageJ hyperstack made by
[ImageStreamCombiner](https://github.com/BioImaging-NKI/ImageStreamCombiner) (one frame per event).

The import page (in the browser, on your own computer): add the zip files, name the channels and give them their
role, set the pixel size and the analysis settings, run all datasets and see an overview of the results.
Every dataset gets a `dashboard.html` with its results (open it from the overview, or from its folder): the event
images with the cells outlined, the classification, the interface measurements and the membrane enrichment, all
linked.

```powershell
uv run imagestream-cli gui
```

The command line:

```powershell
uv run imagestream-cli run-project D:/results/my-project   # the project saved by the page
uv run imagestream-cli example-settings > settings.toml     # or one dataset with a settings file
uv run imagestream-cli run settings.toml
```

## Development

```powershell
uv sync
uv run pytest           # unit tests; the comparison with the Fiji reference needs local test data (data/reference_test.toml)
uv run ruff check
uv run mypy             # strict
uv run bumpver update --patch   # new version: pyproject.toml + version.py, commit and tag
```

## License

GPL-3.0-or-later. Bram van den Broek and Rolf Harkes, The Netherlands Cancer Institute (NKI), Amsterdam.
