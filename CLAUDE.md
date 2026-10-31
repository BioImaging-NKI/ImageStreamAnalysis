# CLAUDE.md

Notes for Claude Code when working in this repository. Last updated 2026-10-22.

## This project: the Python port (started 2026-10-22)

This folder is a **Python port of Fiji macro v1.8**, modelled on FAST-HIPPOS-Python (see "Design inspiration" below).

- Not a clone of upstream: the upstream macros stay in https://github.com/BioImaging-NKI/ImageStreamAnalysis.
- **No remote yet: do not push** until the user asks.

### Confidential research information: never commit it

This repository may become public. **Never put research information in a tracked file, a commit message or a test**:
- the location of the test data, and sample, donor, cell-line, marker, antibody or channel names of real data;
- numbers measured on real data (counts, intensities, enrichment values, IoUs, normalization factors, ...);
- local paths of this machine or of colleagues.

Keep these in `CLAUDE.local.md` (gitignored, read by Claude Code) and in `data/` (gitignored). Tests on real data
read their paths and settings from a gitignored file (`data/reference_test.toml`) and compare with the reference
files instead of hard-coding numbers. Before every commit, check the staged diff for such information. Synthetic
data and generic statements ("matches the Fiji reference to 3 decimals") are fine.

A local guard backs this up:
- `.git/hooks/pre-commit` and `.git/hooks/commit-msg` refuse a commit whose staged files, file names or message contain a term from `data/forbidden_terms.txt`. Both the hooks and the term list are local and not committed.
- Add new sample, marker or path names to that list when new data arrive.
- Never bypass it with `--no-verify`.
- After a fresh clone, the hooks must be installed again; their source is described in `CLAUDE.local.md`.
- The dashboard design is **deliberately postponed**: settle the pipeline first.

### Layout

- `src/imagestream_analysis/` (uv, hatchling, Python 3.12+; `.python-version` is 3.13). Modules:
  - `config.py`: Settings dataclasses ↔ TOML. Relative paths resolve against the settings file, and unknown keys are an error. Channels are numbers or names; `channels.immune` and `channels.functional` are lists; `ChannelSettings.migrate` reads the old `t_cell`.
  - `io.py`: `read_stack` reads only the first N events; also the labelmap montage↔stack and the TSV writer (3 decimals, `NaN`).
  - `imagej.py`: exact re-implementations of the ImageJ/CLIJ2 operations (Otsu on float, RankFilters variance, cross erode/dilate, close index gaps).
  - `segmentation.py`: marker normalization, Cellpose input montage, Cellpose subprocess, label clean-up.
  - `_cellpose_runner.py`: a standalone script that must **not** import the package; it handles Cellpose 3 and 4.
  - `measure.py`: the first stage: the per-cell measurements (normalized markers, every loaded channel minus its background, area in µm²) and cells cut by the edge of their event image.
  - `classification.py`: k-means on the measured cells, tumor/immune-cell assignment, class map.
  - `interfaces.py`: membrane and interface masks, the per-event table, doublets.
  - `membrane.py`: the optional membrane analysis (`membrane.analysis`): outline, straighten, profiles, unseeded bootstrap control, aligned average, histograms.
  - `pipeline.py`: `load(settings)` reads an IDEAS zip or a Combiner `.tif`; `run(settings)` writes the output folder.
  - `ideas.py`: reads IDEAS zips directly (replaces ImageStreamCombiner).
    - Events are ordered by object number, and the numbers are kept as an `Object ID` column.
    - Only the needed channels are loaded.
    - It also makes the channel thumbnails.
    - The OME metadata (operator, source path) is never read or shown.
  - `project.py`: `project.toml` (the panel: roles by marker name; zips with their channel names; datasets with experiment labels; shared settings). Old projects with roles by number per zip are converted on loading.
    - `run_project` analyses every dataset into `<project>/<experiment>/<dataset>/` and writes `summary.json` and `overview.tsv`.
  - `dashboard.py` + `templates/dashboard.html`: the results page of each dataset (`dashboard.html` + `dashboard_files/*.webp`).
    - Written by `pipeline.save`.
    - Self-contained and works from file://. Event images are drawn onto canvases but never read back (tainted canvases are fine).
    - Linked hover and selection; `#cell=N` / `#event=N` in the address selects on opening.
  - `gui.py` + `templates/gui.html`: the import and run page (`imagestream-cli gui` / `imagestream`).
    - A stdlib HTTP server on 127.0.0.1 with a token; the page is plain HTML/JS with no external resources.
    - Python, not Node: the server reads zips and starts the analysis.
  - `cli.py`: the `imagestream-cli` entry point.
- `reference/fiji_v1.8/`: the macro, copied from upstream commit `cffe8fb`. Only v1.8 matters.
- `docs/SPEC.md`: the behaviour of the macro with line numbers, port status per step (✅/🟡/⬜), and the issue list I1–I14. **Update it with every port change**, and record deviations there with their reason.
- `tests/`: `conftest.py` (synthetic touching-cell events), `test_imagej.py`, `test_pipeline.py`, and `test_reference.py`.
  - `test_reference.py` compares with the Fiji output of the (confidential) test data; it reads `data/reference_test.toml` (or the file in `ISA_REFERENCE`) and is skipped without it.
- `data/`, `output/` and `CLAUDE.local.md` are gitignored.

### Commands

- `uv sync`
- `uv run pytest`: about 30 tests in under 10 s, including the reference comparison when its data are available.
- `uv run ruff check`: line length 120, rules E, F, UP, B, SIM, I (taken from `rharkes/lif2ometiff`).
- `uv run mypy`: **strict**; it must stay clean. `sklearn`, `cellpose` and `torch` are ignored as untyped; pandas-stubs and scipy-stubs are installed.
  - Pattern for returned numpy expressions: assign to a typed local first (`out: np.ndarray = ...; return out`).
- `uv run bumpver update --patch|--minor|--major`:
  - It updates `pyproject.toml` and `src/imagestream_analysis/version.py`, then commits and tags.
  - `push = false` until there is a GitHub remote.
  - Versions are plain `MAJOR.MINOR.PATCH`, with no `.dev` suffix.
  - It needs a clean working tree with at least one commit.
- `uv run imagestream-cli gui` opens the import page; `run-project <folder>`, `run <settings.toml>` and `example-settings` are the other commands.
- To check the page visually: start `imagestream-cli gui --no-browser --port 8765`, then take a screenshot with headless Edge (`msedge --headless=new --screenshot=... http://127.0.0.1:8765/`). The page loads the last project from `~/.imagestream_analysis/gui_prefs.json`.

### Validation status (2026-10-22)

- Against the Fiji reference output of the test data (details and numbers in `CLAUDE.local.md`):
  - Cells and event tables: match to 3 decimals, 100 % class agreement.
  - Membrane analysis: validated on distributions.
- Lessons learned:
  - **ImageJ Otsu on 32-bit:** byte = `round((v-min)*255/(max-min))`, and lower = `min + (k+1)/255*(max-min)`. Other variants miss by 0.3–1.5 %.
  - **Cell order within an event is nondeterministic in Fiji,** so compare cells matched on (event, area). The reference table has no event column; derive it from the per-event counts.

Done (2026-10-22): reading IDEAS zips, the import page, and the dashboard (v1, to be refined with the user).

### Decisions by the user (2026-10-22)

- **The port need not replicate Fiji exactly.** The reference output is a sanity check, not a target. Improving on the macro is fine, as long as each deviation is recorded in SPEC.md with its reason and its effect on the test data.

- **Segmentation (I2): always a montage.** Single event images don't look like the microscopy images Cellpose is made for. Do **not** propose per-event segmentation again; it was tried and was slower and worse.
  - The montage need not be a strip: the default is an about-square grid (`segmentation.montage_columns = 0`), and `montage_columns = n_events` gives the macro's single row.
  - Cells crossing a tile border are split per event (I15).
- **Tumor cluster (I4):** keep the macro's *sum* of ratio/area.
- **Membrane analysis, decided 2026-10-22 (next port step):**
  - **I8:** one membrane per cell, not per blob of touching same-class cells.
  - **I9:** mark the end of each profile explicitly; true zeros are data.
  - **I10:** the control is **not** seeded. Cells without an interface get `control_repeats` (10) random windows, as long as the median interface length; all windows are pooled for the control distribution. Per-distribution statistics wait for the dashboard.
  - The membrane analysis is optional (`membrane.analysis`, off by default): not every ImageStream analysis needs it.
  - **I11:** the average profile only uses cells with an interface.
    - Verified on the test data that the macro's average includes the rings without an interface, which roughly halves the peak (numbers in `CLAUDE.local.md`).
  - **Validation:** compare distributions with Fiji (`stdDevs above mean`, histograms, the shape of the average profile), not cell by cell.
- **IDEAS import page (done 2026-10-22):**
  - It is a separate local page.
  - It takes several zips from several experiments, with channel names and roles and the pixel size per zip ("Copy setup to all").
  - It shows an overview across experiments.
  - A Node.js server was considered and not used: the server must read zips and run Python anyway.
- **Dashboard (v1 done 2026-10-22):** one `dashboard.html` per dataset, opened from the overview of the import page.
  - Colours come from the dataviz skill's validated palette: T cells blue, tumor orange, interface aqua.
  - The design follows the FAST-HIPPOS mapping table below. Not yet built: the cross-dataset overview page, and live re-classification.
- **I6 was a false alarm.** Interface XOR membrane equals membrane AND NOT interface, because the interface lies inside the membrane band (same k, same cross element). Check claims like this on the data before raising them.

### Decisions by the user (2026-10-31): merging `ClassifyCellTypes_2_8.ijm`

The later macro (several immune subtypes, cell interactions, no membrane analysis) is merged into this tool, so that
v1.8 is the special case "one immune marker, no functional markers". Not built yet.

- **Where 1.8 and 2.8 differ, pick what is best.** Neither macro has to be mimicked exactly; record the choice and its
  reason in SPEC.md.
- **Normalize before clustering.** 2.8 computes normalized functional-marker values but adds the raw ones (line 650):
  a programming mistake.
- **The "crap" class gets fixed and renamed.** It stands for detections that are not real cells: positive in too
  many marker channels. (In 2.8 it is switched off and cannot work.)
- **Subtype thresholds:** MaxEntropy was chosen because it worked well. It may be replaced by something similar that
  exists in Python, if porting it is a lot of work.
- **Membrane enrichment between two cell types of choice**, and ideally chosen in the dashboard. With several cell
  types the number of possible interfaces grows, which is why the membrane part of 2.8 stayed unfinished.
- **Pooling** of the classification across datasets: left to Claude's judgement.
- **Cellpose montage: a grid** (as already decided for 1.8).
- Storing and reusing settings, and several inputs per run, are part of the design.
- The IDEAS exports are compensated (as far as the user knows), and dead cells are already removed in IDEAS. A
  viability-marker role may come later: keep room for it.
- Flag (and optionally exclude) cells cut by the edge of the event image.
- At most 8 cell types (the 8 categorical palette slots); artefacts in grey.
- The poster is left as it is for now.
- A user manual (MkDocs, as FAST-HIPPOS) will be started once the new settings have settled.
- **Classification stays in two steps** (lineage, then immune subtypes), as in 2.8; see SPEC.md.

### Git

The history of the first machine is not here (2026-10-31): the repository was started again with a snapshot commit
of the same files. To join the two later: fetch the old repository as a remote and
`git rebase --onto <old main> <snapshot commit> main` (the snapshot's tree equals the old tip, plus `.gitignore` and
`.python-version`).

## General environment

- Windows 10, PowerShell 5.1 and Git Bash. `python` on PATH is only the Microsoft Store stub, so always use `uv run`. `pixi` and Node are also installed.
- The IDE's type checker uses a different interpreter, so its "Cannot find module numpy" diagnostics are noise.
- Use metric units by default (organization policy).

## Pipeline overview

The data come from the **Amnis/Cytek ImageStream**, an imaging flow cytometer. They do **not** come from the ImageXpress plate imager; `BioImaging-NKI/ImageXpress_converter` is unrelated to this project.

1. Export single-cell images from IDEAS as `<objectID>_Ch<N>.ome.tif`, one subfolder per dataset, and zip them.
2. **ImageStreamCombiner** (https://github.com/BioImaging-NKI/ImageStreamCombiner) turns the zip into one ImageJ hyperstack `.tif` per dataset.
3. **ImageStreamAnalysis** (upstream repo) runs Fiji macros that segment and classify cells and measure interactions, then a notebook summarizes the results. This folder ports that step to Python.

## ImageStreamAnalysis (upstream repo with the Fiji macros)

Authors: Bram van den Broek (b.vd.broek@nki.nl) and Rolf Harkes. The README says the repo is "mostly a placeholder"; there are no run instructions, tests or CI.

### Files

- `Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm` (about 1280 lines): **the version to use.** This is the two-cell-type analysis, added 2026-10-22 ("Speed improvement and bugfixes"). Work only on 1.8.
- `Analyze_..._1.5.ijm`: an obsolete older version that is still in the repo. Do not use it, edit it or take code from it. Versions are kept as separate files, not in git history.
- `ClassifyCellTypes_2_8.ijm` (about 1680 lines): a more general classifier for several immune subtypes. It does not do membrane enrichment analysis. The README says the two macros are meant to be merged eventually.
- `Summarize_Interactions.ipynb`: pandas notebook that cross-tabulates cell-type interaction pairs and draws a stacked bar chart.

### Input format (from ImageStreamCombiner)

- One ImageJ hyperstack per dataset with axes TCYX: one **frame per ImageStream event** (cell image).
- Each channel's slice label holds the channel name.
- The `Medians` property is written as `"\n" + "m1,m2,..."`, so the macros read the comma-separated per-channel **background medians** from line 2 of `getImageInfo()`. They are used for display ranges and background subtraction, and padded pixels are replaced with them.
- Macro channel numbers are positions among the channels **named** in the Combiner, not the original ImageStream `Ch` numbers.
- Default channels: 1 brightfield, 2 marker of interest (MOI), 3 tumor marker, 4 T-cell marker.

### Analyze_..._1.8 pipeline

1. Script parameters (`#@`) are saved to `settings.txt` in `<output>/<input name>/`.
2. Display LUTs are set, and a montage of all frames is built.
3. **Segmentation** with Cellpose (default `cyto3`, diameter 30 px, flow threshold 0.5, `--use_gpu`) through the PT-BIOP Fiji wrapper. Input is the brightfield variance image plus the summed nuclear markers (if `markerLoc = nuclear`). A saved labelmap can be loaded instead of running Cellpose.
4. Labels are eroded on the GPU (CLIJ2), then filtered to 100–10000 px.
5. **Classification**: per-cell marker intensities are z-scored, then k-means with 2 clusters and fixed seed 48. The cluster with the lowest marker ratio is called tumor.
6. **Interfaces**: membrane rings (`membraneThickness` px), the contact zones between touching cells, and the non-contact parts are built.
7. **Outputs** (`.tsv` and `.png`): tables `All cell combinations`, `Cell classification`, `Single Tumor Single T cell` (doublets) and `Membrane intensity analysis`, plus histograms.
8. **Membrane enrichment**: membranes are straightened, and MOI intensity is normalized to the non-interface part (mean 0, SD 1). A control interface at a random position is used for comparison. Profiles are aligned so the interface is centered at position 20, then averaged.
- The Cellpose `env_path`/`env_type` are read from the BIOP wrapper's Fiji prefs rather than from the dialog, so set them once in the BIOP Cellpose command. An old `Cellpose Advanced` wrapper makes the macro exit.

### ClassifyCellTypes_2_8 pipeline

- Takes several input stacks (`File[]`) and a `baseName` (default `Panel_1`).
- Uses phenotypic and functional immune marker channels, given as comma-separated lists.
- Tumor vs immune is split with k-means (seed 48). Immune markers are then thresholded positive/negative (automatic or manual threshold) to assign an `Immune subtype`.
- It records which cells touch through a region adjacency graph, stored in the `Interaction with (IDs)` column as 1-based IDs.
- Outputs: `<baseName>_Cell statistics.tsv`, `<baseName>_Image statistics.tsv`, `<baseName>_multiplets.tsv`, `<baseName>_labelmap_cells.zip` and `<baseName>_celltype_clustering.tif`.
- Known caveat (from a code comment): background subtraction uses only the medians of the **last** input image.

### Summarize_Interactions.ipynb

- Reads a hard-coded `Measurements.tsv` (ISO-8859-1, tab-separated). This name does not match what the macro saves (`<baseName>_Cell statistics.tsv`), so the file must be renamed or the path edited.
- Columns used: `Class` (Immune/Tumor), `Immune subtype`, `Interaction with (IDs)`. Subtypes `Unknown` and `All negative` are mapped to `Other`.
- The interaction IDs index rows directly (`df.iloc[id-1]`), so the TSV must hold one row per cell, in ID order.

### Fiji dependencies

Update sites: CLIJ, CLIJ2, IJPB-plugins (MorphoLibJ) and PTBIOP (Cellpose wrapper, "Label image to composite ROIs"). The k-means command also needs its plugin. A working Cellpose Python environment is required, and a GPU is expected.

### Conventions in the macros

- Use ImageJ macro language with `#@` script parameters.
- Constants are UPPER_CASE at the top of the file; variables are camelCase.
- Heavy work goes through `Ext.CLIJ2_*` on the GPU, with `setBatchMode` used for speed.
- Every run saves its settings with the results.
- Utility functions are at the bottom of the file. In 1.8, several helpers (`firstIndexOfArray`, `sumArray`, ...) are defined **twice**.

## Design inspiration: FAST-HIPPOS-Python

https://github.com/Jalink-lab/FAST-HIPPOS-Python (GPL-3.0, Bram van den Broek, actively developed in 2026-10). It is a Python port of a Fiji macro by the same group, so it is the template for how this project should look.

### Architecture worth copying

- **Packaging:** uv with `uv.lock`, src layout (`src/fast_hippos`), hatchling and Python 3.12 or later.
  - Two entry points: `fast-hippos-cli` (console) and `fast-hippos` (gui-script, no console).
  - Dev group: pytest, ruff (line length 110) and mini-racer.
  - The manual is MkDocs Material in `manual/` and is published to GitHub Pages by a workflow.
- **Settings:** dataclasses in `config.py`, loaded from TOML; they replace the Fiji `#@` dialog. Channels are 1-based, as in Fiji. A run saves its full settings TOML next to the output, which fixes the gap in `settings.txt` here.
- **Reference implementation:** the Fiji macro is kept in `reference/fiji_v0.9.5/`. `docs/SPEC.md` describes its behaviour with macro line numbers, together with deliberate deviations and the decision date for each. It also validates the port against Fiji output.
- **Cellpose:** `_cellpose_runner.py` is a standalone script, called as `python runner.py in.npy out.npy '{json}'`. It can therefore use the Python of **another environment** (Cellpose 3 or 4), and it reports the version, device and timing.
- **GUI:** a local web app (`gui.py` with `templates/gui.html`).
  - It uses a stdlib `ThreadingHTTPServer` on 127.0.0.1 with a per-session token.
  - The settings form is generated from the dataclasses, with a `FIELDS` dict for labels, help texts, basic/advanced levels and show-conditions.
  - Runs start as a subprocess, and its log is streamed to the page.
  - State is kept in `~/.fast_hippos/`. A second start reuses the running instance, and the GUI stops after 15 minutes idle.
- **Outputs are self-contained HTML:**
  - Each image gets a `dashboard.html`; each run gets an `index.html` overview.
  - `dashboard.py` assembles them from `templates/dashboard/*.js|css` and `templates/common/`, and inlines the data as JSON (float arrays as base64, images as PNG data URLs).
  - They work from `file://` with no server, so they can be sent to a colleague.
- **Live re-analysis in the browser:** `templates/screening.js` recomputes the hit criteria in JS. `tests/test_dashboard_js.py` runs it in mini-racer and checks parity with the Python code.
  - Settings changed in the dashboard can be exported, then **re-applied** to the run with `reapply` without segmenting again.
  - An image can have its own settings that override the run settings.
- **Tests:** synthetic data generator (`synthetic.py`), pipeline tests, GUI tests, JS parity and syntax tests, and manual checks.

### Dashboard UX

- **Top bar:** name, counts and event chips, frame player, LUT and range, view options, box/lasso select, a link up to the overview, and a light/dark toggle.
- **Panels:** Image (zoom/pan, outlines of all cells, hits or the selection; a scale bar), Traces (lines, or median with bands, or density), Kymograph (sortable), Distribution, Scatter, Cell card, Table, and a foldable Screening pane with live criteria.
- **Everything is linked:** hovering a cell anywhere highlights it in every view, and selections (click, Ctrl-click, Shift-drag box or lasso) are shared. Keyboard shortcuts: O, Z, F, 1, L and Esc.
- **Overview:** one row per image with thumbnail, counts, medians and **QC flags** (warnings only). Clicking a name opens the image's dashboard. Export works per selection.
- **Help:** a `?` on each panel opens the matching part of the manual (`help.json`) in a side panel.

### How the concepts map to ImageStream (proposal, not yet agreed)

| FAST-HIPPOS | ImageStreamAnalysis |
|---|---|
| One large image with cells | A **gallery** of events, one tile per frame. Channel composite, outlines coloured by class, interfaces drawn |
| Trace per cell over time | **Membrane profile** per cell, aligned so the interface is at 20, plus the average |
| Kymograph (cells × time) | `membrane_profiles_normalized` (cells × perimeter position), sortable by enrichment |
| Scatter / distribution | Marker ratio vs area (tumor/T-cell classes); interface vs non-interface MOI; histograms of interface vs control SDs |
| Live screening criteria | Live class threshold or cluster override, and filters such as "single tumor + single T cell" or size gates |
| Run overview with QC | Samples/donors side by side. QC examples: few doublets, unbalanced clusters, no enrichment vs control |
| LIF/TIFF readers | Read the IDEAS **zip directly** (absorbing ImageStreamCombiner) and **keep the object IDs** |

## ImageStreamCombiner (upstream repo)

https://github.com/BioImaging-NKI/ImageStreamCombiner, v1.0.2, GPL-3.0, Python 3.11 or later, last changed 2025-03. Author: Rolf Harkes. Package `combine_imagestream_files` (src layout, built with meson-python).

### How it works

- `src/combine_imagestream_files/imagestreamzip.py`, class `ImageStreamZip`:
  - `loadfile(zip)` scans the zip. Every **first-level subfolder** is a dataset, and only files ending in `.ome.tif` that sit directly in one count. Files in the zip root or deeper are ignored.
  - `writetiffs(folder, pixelsize)` builds the output per dataset:
    - A cell is skipped (with an error logged) if it lacks any named channel.
    - Each cell gets one median per channel, and the dataset median is the median of those.
    - Images are padded to the largest width and height, centered, with the dataset median as fill.
    - The data type comes from the first file read.
    - It writes `<dataset>.tif` with `imagej=True`, axes TCYX, pixel size in µm, channel `Labels`, display `Ranges` (median to max) and the `Medians` property.
- `dataset.py`:
  - `DataSet.addfile` parses `<stem>_Ch<N>` and groups files by stem (the cell).
  - `Channel` has a 1-based index and a name. An empty name, or `--` in the GUI, means the channel is not exported.
- `groupthem.py` holds an older standalone version of the same grouping; it is unused.

### Front ends

- CLI entry point `combine_imagestream_files -i <zip> -p <pixelsize> -l 0|1|2`.
  - With no channel names set, the CLI exports **nothing**: every channel is invalid and it logs "No valid channels". In practice only the GUI is usable.
- `ImageStreamCombiner.py`, a PySide6 GUI:
  - You can drop a zip or a `.toml` file on the window.
  - Channel names can be edited, saved to and loaded from `.toml` (`[channels] 1 = "DAPI" ...`, see `demofiles/Arctic_Conditions.toml`), and copied to all datasets.
  - The pixel size is set in a spin box.
  - Output goes next to the zip.
  - Settings are stored in QSettings under `NKI-AVL/CombineImageStreamFiles`.
  - It writes `debug.log` to the current working directory.
- `ui/form.ui` is converted with `pyside6-uic -o ui/form.py ui/form.ui`. Never edit the generated `form.py`.

### Build and tooling

- Editable install: `pip install meson-python meson ninja`, then `pip install --no-build-isolation --editable .[dev]`.
- Executable: `pyi-makespec` and then `pyinstaller`.
- CI: black, plus mypy strict on Python 3.10–3.12. Note that the code uses `Path | str`, which needs Python 3.10 or later, while `requires-python` says 3.11 or later. There are no tests.
- `demofiles/Demo.zip` holds the dataset `Arctic_Conditions` with 12 channels.

### Gotchas that affect the analysis

- **Frame order is the zip listing order, not object-ID order.** The cell stems (IDEAS object IDs) are **not written** to the output, so a frame cannot be traced back to its IDEAS object.
- **Channel numbering changes:** only named channels are exported, in index order. With `Arctic_Conditions.toml` the stack has 5 channels: DAPI, PanCK, CD3e, CD4, CD8, with **no brightfield**. The macro defaults (BF = 1) assume a brightfield channel was named.
- The pixel size defaults to 1.0 in the CLI. Set it to the real value; nominal ImageStream values are about 1 µm at 20×, 0.5 µm at 40× and 0.33 µm at 60×, but verify them for your instrument.
- **Memory:** the whole dataset is held in RAM as a numpy array, and every file is read twice.
- **GUI:** writing runs on the UI thread, so the window freezes while it works.
