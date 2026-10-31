# ImageStream Analysis: behaviour specification for the Python port

Reference: Fiji macro **v1.8**,
`reference/fiji_v1.8/Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm`
(BioImaging-NKI/ImageStreamAnalysis, commit `cffe8fb`, 2026-10-22). Line numbers below refer to that file.
Older versions of the macro (1.5) are not used.

Validation data: a confidential test data set with the output of a v1.8 run (not in this repository; see
`tests/test_reference.py`). Numbers measured on it are kept out of this file.

Status per step: ✅ ported and validated against the reference · 🟡 ported, not validated · ⬜ not yet ported.

---

## 1. Input (lines 132-150) ✅

- ImageJ hyperstack from ImageStreamCombiner, axes TCYX, one frame per ImageStream event (`io.read_stack`).
- Channel names: the slice labels. Pixel size: the TIFF resolution.
- Background per channel: the `Medians` property (`"\nm1,m2,..."`); the macro reads line 2 of `getImageInfo()`.
  If missing, the macro uses an empty array (zeros); the port uses 0 and logs a warning.
- `limitNrCells` (`input.max_events`): only the first N events. The port reads only those pages from disk.
- The macro converts the whole stack to 32-bit (line 157); the port converts per channel when it is used.

**IDEAS zip, read directly (new; replaces ImageStreamCombiner, `ideas.py`)** ✅

- An IDEAS export: a zip with one folder per dataset, one file per object and channel, `<object>_Ch<N>.ome.tif`
  (uint16, sizes differ per object). The files hold no pixel size and no channel names (`input.pixel_size`,
  `input.channel_names`, from the user). Their OME metadata (operator, source path) is not read.
- Only the channels with a role plus `input.channels` are loaded; the channel settings then use the IDEAS numbers.
- Objects without all loaded channels are left out (logged). Events are ordered by **object number** and the numbers
  are kept: an `Object ID` column is added to the cell, event, doublet and membrane tables.
  (The Combiner orders by zip listing and drops the numbers, so `max_events` selects different events.)
- Background per channel and padding as the Combiner: the median over the objects of their median; every image is
  centred in the size of the largest one, padded with that background. Verified on the test data: the same
  backgrounds, size and pixels as the Combiner stack.

## 2. Marker normalization (`normalizeImage`, lines 1145-1160) ✅

For the tumor and T-cell marker channels (and in membrane mode also the brightfield variance):

1. Per event: `setAutoThreshold("Otsu dark")` on the 32-bit frame, `getThreshold(lower, upper)`.
   ImageJ does this on an 8-bit version: display range reset to min-max, byte = `round((v-min)·255/(max-min))`,
   Otsu (AutoThresholder) on the 256-bin histogram gives bin k, lower = `min + (k+1)/255·(max-min)`.
   **Verified**: this variant reproduces Fiji's normalization factors on the test data; flooring or 256-scaling
   miss them by 0.3-1.5 %.
2. Factor = `medianOfArray(thresholds)` = element `floor(n/2)` of the sorted values (upper median).
3. The channel is divided by the factor. `Enhance Contrast` only changes the display.

## 3. Segmentation (lines 199-279) 🟡

Cellpose input, all events side by side in **one row** (Make Montage, rows=1):

| `markerLoc` | Cellpose image | Cellpose 3 channels |
|---|---|---|
| nuclear | C1 = Variance r=2 of brightfield (not normalized), C2 = tumor + T-cell marker (normalized) | ch1=1, ch2=2 |
| membrane | tumor + T-cell marker + brightfield variance, each normalized | ch1=1, ch2=0 |

- Variance filter: ImageJ RankFilters, circular kernel (radius 2: 5×5 without corners, 21 px), population variance,
  edges padded with the nearest pixels (`imagej.variance_filter`).
- Cellpose: model `cyto3`, diameter (`cellDiameter`), flow threshold 0.5 (constant), cellprob 0, `--use_gpu`.
- Or `open(labelmapFile)` (`segmentation.method = "labelmap"`). The macro saves the Cellpose montage as
  `labelmap_cells.zip` only when Cellpose ran (line 284); the port saves `labelmap_cells.tif`. Both are accepted.

**Port:** Cellpose runs in a subprocess (`_cellpose_runner.py`) with the Python of any environment
(`segmentation.cellpose_python`), Cellpose 3 or 4. With Cellpose 4 the model is `cpsam`, which takes the channels as
they are (no `channels` argument); the diameter is passed on (cpsam rescales by 30/diameter; 0 = no rescaling).
**Grid montage (default, decision on I2):** the events are put in an about-square grid, row by row
(`io.grid_columns`), empty tiles at the end are 0. The labelmap is saved as the montage (`labelmap_cells.tif`);
`read_labelmap` recognizes a row or a grid.

Validated with Cellpose 3 `cyto3` on the test data (2026-10-22): as one row (`montage_columns` = number of events,
as the macro) it reproduces the Fiji labels almost exactly; the grid is about 30 % faster and gives slightly
different masks (median IoU with the Fiji cells about 0.97). Per-event segmentation was tried and dropped: slower
and further from the reference. `cpsam` (Cellpose 4) runs but segments differently: a model choice, not a port
difference.

## 4. Label clean-up (lines 281-317) ✅

1. `closeIndexGapsInLabelMap`: labels renumbered 1..N in their order.
2. If `erosionSteps > 0`: binary erosion of the foreground (all labels together), `erosionSteps` times
   `erodeSphereSliceBySlice` (4-connected cross; outside the image = nearest edge pixel), multiplied with the labels.
   Then `excludeLabelsOutsideSizeRange(minCellSize=100, maxCellSize=10000)` (kept if min ≤ size ≤ max) and gaps
   closed again.
3. Cells are made into ROIs ("Label image to composite ROIs", BIOP) and measured with the ROI Manager.

Cell order: by event, then label. In Fiji the order within an event is **not deterministic** (issue I1).

**Cells cut by the edge of their event image (new, 2026-10-31)** ✅

Each event's own image is smaller than the frame when it was padded (by ImageStreamCombiner or the IDEAS reader) to
the size of the largest event. Its extent (`EventStack.image_extents`) is known when read from an IDEAS zip; in a
Combiner stack it is found as the outer rows and columns that hold exactly the background value in every channel
(on the test data the same sizes as the zip gives). A cell within `erosion` pixels of the edge of its event's image
(it touched the edge before the cells were shrunk) gets `Edge` = 1. `segmentation.exclude_edge_cells` (default off)
leaves them out, before the classification, with the gaps in the label numbers closed.

## 5. Per-cell measurements and classification (lines 329-507) ✅

The port runs in two stages (`pipeline.measure`, then the rest of `pipeline.analyze`): segmentation, clean-up and
per-cell measurements (`measure.measure_cells`) first, the classification and everything after it on the measured
table. New measurements, as in ClassifyCellTypes 2.8: `Area (µm²)`, `Edge`, and `Mean ch<N> (<name>)` for every
loaded channel: the mean of the stored values minus the channel's background (median), not normalized.

- `Area` = pixel count; `T-cell intensity`, `Tumor intensity` = mean of the **normalized** marker channels.
- `Marker ratio` = T-cell / tumor intensity; `Marker ratio over Area` = ratio / area.
- k-means (ImageJ "k-means Clustering", 2 clusters, tolerance 1e-4, seed 48) on three features per cell:
  log(T-cell intensity), log(tumor intensity), area, each z-scored (mean 0, sample sd).
  Port: scikit-learn `KMeans(n_init=10, random_state=48)` (issue I14).
- Tumor cells = the cluster with the lower **sum** of `Marker ratio over Area` (issue I4).
- Class map: T cells → 1, tumor cells → 2 (`changeValues`). A third class (`class2`, green) cannot occur with 2
  clusters.
- `Class` column: in the macro the arbitrary cluster number; in the port **1 = T cell, 0 = tumor cell**
  (on the test data this equals the reference). The port adds `Frame`, `Label` and `Cell type` columns.
- Plot `cell Classes.png`: area vs marker ratio, log-log, T cells blue, tumor cells orange. ⬜ (dashboard)

## 6. Membranes and interfaces (lines 509-533) ✅

With k = `membraneThickness`, per class mask (T cells / tumor cells, all cells of a class together):

- membrane = mask AND NOT erode(mask, k) (4-connected, k times)
- T-cell interface = T-cell mask AND dilate(tumor mask, k); tumor interface likewise
- "membrane not in interface" = interface **XOR** membrane. Because the interface lies inside the membrane band
  (see I6), this equals membrane AND NOT interface.

The overlays on the montage (lines 555-631) are display only. ⬜ (dashboard)

## 7. Per-event measurements (lines 633-698) ✅

`All cell combinations`: one row per event: `Frame` (1-based), `T cells`, `Tumor cells`, and the mean of the
**raw** marker-of-interest channel (32-bit, no background subtraction) within each of the four masks of that event
(NaN when empty), plus `T cell membrane ratio` = interface / not-in-interface and `Tumor membrane ratio`.

`Single Tumor Single T cell`: the rows with 1 T cell, 1 tumor cell and a non-NaN T-cell interface.

Plots of interface vs non-interface MOI for the cell of interest (all events, doublets). ⬜ (dashboard)

## 8. Straightened membrane profiles (lines 700-827) ✅ (distributions; optional: `membrane.analysis`)

For the cell of interest (`CellOfInterest`; both `markerLoc` branches are identical):

1. Particles of the membrane mask (Analyze Particles, `add stack`): one ROI per connected membrane ring, i.e. per
   blob of touching cells of that class (issue I8). On the test data this gives about 10 % fewer rings than cells.
2. MOI minus its median background (line 742).
3. `create_straightened_membrane_MIP_profiles`: per ROI, Area to Line (traced outline), `Straighten` with line width
   `membraneThicknessForStraightening` (5), maximum projection across the width → one profile per ring;
   `traceLengths` = straightened length. Profiles left-aligned in an image (rows = rings), padded with 0.
4. Zeros → NaN (lines 748-750; also true zeros of the data, issue I9).
5. The interface mask is straightened in the same way; > 0 = interface, the rest = non-interface.
6. Per ring: interface mean, non-interface mean and sd, ratio, `stdDevs above mean` = (interface − non-interface
   mean) / non-interface sd. Without interface: a **control** interface of 25 px at a random position
   (`round(random·length − 1)`, unseeded, issue I10) → `control stdDevs above mean`.
7. The profiles are normalized to the non-interface mean and sd → `Membrane intensity analysis.tsv`.

## 9. Histograms (lines 829-891) ✅

`stdDevs above mean` (interface) and `control stdDevs above mean`: `round((maxHistValue+5)/binWidth)` bins over
[−5, maxHistValue] (100 bins with the defaults) → `Histogram_*_stdDevs_above_mean.tsv/.png`.

## 10. Aligned profiles (lines 898-996) ✅

- Trace lengths re-measured as the first NaN of each row.
- Shift per row = 20 − (first interface pixel, or the start of an interface that wraps around x = 0) −
  interface length / 2; all four profile images rotated by it (`Array.rotate`). Rows without interface are rotated
  too (issue I11).
- `membrane_profiles_normalized.tif`; average profile over all rows up to round(mean + sd of the trace lengths) →
  `Average normalized membrane profile values - interfaces centered at 20.tsv` and its plot.

## 8-10 in the port (`membrane.py`)

Optional (`membrane.analysis`, default off: not every ImageStream analysis needs it).

- **One membrane per cell of interest (I8)**: the outline of the cell (marching squares at 0.5 on its own mask, i.e.
  along the pixel edges; smoothed along the line with a Gaussian, σ = 2 samples, in place of Fiji's spline fit;
  resampled at 1 px), straightened with `membrane.straightening_width` samples across it (bilinear), maximum
  across → profile. The cell's interface (the interface mask of its class, within this cell) is straightened along
  the same line; > 0 = interface. Background: the MOI median of the stack. The profile lengths are kept (I9).
- Per cell: interface mean, non-interface mean and sd (sample sd), ratio, `stdDevs above mean`, and the new
  `interface length` column.
- **Control (I10)**: for cells without interface, `membrane.control_repeats` (10) windows at random positions,
  unseeded, of the **median interface length** of the cells with an interface (the macro: one window of 25 px).
  The table gets the mean of a cell's windows; all windows are pooled in
  `Control windows stdDevs above mean.tsv` and in the control histogram.
- **Alignment**: the centre of the *longest* interface run (circular) is put at x = 20 (the macro: first interface
  pixel + total interface length / 2). Cells without interface are not shifted and **not averaged (I11)**.
- Average profile over the cells with an interface, up to round(mean + sd) of their lengths, NaN-aware.
- Histograms: as the macro, but the first column is called `bin_start` (the macro's `bin_center` holds the bin starts).
- `membrane_profiles_normalized.tif`: one row per row of the membrane table, NaN beyond the length. Plots ⬜ (dashboard).

Validation on the test data (2026-10-22), on distributions: the number of cells with an interface and the medians
of the membrane length, the interface and non-interface means and their ratio are within 1-10 % of Fiji's; the
quartiles of `stdDevs above mean` within 0.5 SD; the pooled control windows have the quartiles of Fiji's single
controls within 0.25 SD; the average profile of the cells with an interface correlates > 0.99 with Fiji's profiles
of the rings with an interface. (These are the checks of `tests/test_reference.py`.)

## Projects and the import page (new) ✅

- `project.py`: a project = several IDEAS zips (of one or more experiments) with shared analysis settings, saved as
  `project.toml` in the project folder. Per zip: path, pixel size, channels to load and their names, the roles
  (brightfield, marker of interest, tumor, T cell), and per dataset an experiment label and whether to include it.
  `imagestream-cli run-project` analyses every included dataset into `<project>/<experiment>/<dataset>/`; a failing
  dataset is logged (and recorded in its `summary.json`) without stopping the others. `overview.tsv` holds one row
  per dataset (counts, doublets, median membrane ratio of the doublets, and the membrane results if computed).
- `gui.py` + `templates/gui.html` (`imagestream-cli gui`, or `imagestream`): a local page to add zips, name the
  channels with example images of each, give them a role, set the pixel size and the analysis settings, run with the
  log streamed, and see the overview. Server: stdlib `ThreadingHTTPServer` on 127.0.0.1 with a per-session token;
  stops after 15 min without requests; a second start reuses the running one.

## The dashboard (new; replaces the macro's plots and overlays) ✅

`dashboard.py` + `templates/dashboard.html`: every analysed dataset gets `dashboard.html` (results embedded as JSON)
and `dashboard_files/` (event images: per channel, 8-bit, display range background (median) to the 99.8th
percentile, brightfield 0.2-99.8th percentile; sheets of 10 x 10 events, WebP). It opens from disk (file://) or
from the import page (served from the project folder, token-protected). No external resources.

- Gallery of the events (only the visible tiles are drawn): channels composed in the page (colour per channel,
  on/off, brightness; brightfield starts at half brightness), cell outlines in the class colour, interfaces.
- Cell classes: area against marker ratio, log-log (the macro's `cell Classes.png`).
- Per event: MOI at the interface against the rest of the membrane of the cell of interest (the macro's
  `Channel N in ... membrane plot`), with the identity line.
- Membrane analysis (if computed): histograms of interface and control as fractions, the average aligned profile,
  and a heatmap of all profiles (diverging, ±6 SD).
- Cell card, cell and event tables (sortable). Hover and selection are linked in all views; one filter row (events:
  all / both cell types / one tumor + one T cell / with an interface; sort) scopes everything. `#cell=N` or
  `#event=N` in the address selects on opening.
- Colours: the validated reference palette (dataviz skill): T cells slot 1 blue, tumor cells slot 2 orange,
  interface slot 3 aqua, in light and dark mode; aqua is below 3:1 on the light surface, so it is always named in
  a legend and the values are in the tables.
- Cost on the full test dataset: about 6 s; outlines are traced on each cell's bounding box.

## Merging ClassifyCellTypes 2.8 (planned, decided 2026-10-31)

The later macro `ClassifyCellTypes_2_8.ijm` (same upstream repository) classifies several immune subtypes and
records which cells touch; it has no membrane analysis. It is merged into this port so that v1.8 is the special case
"one immune marker, no functional markers". Neither macro is mimicked exactly where they differ.

**Classification in two steps, as 2.8 (decision 2026-10-31).** First the lineage: tumor against immune cells
(k-means, 2 clusters, as now); then the immune subtypes, by a threshold per functional marker, computed on the immune
cells only (2.8 thresholds on all cells). Not one clustering of all cells into k types, because:
- k-means favours clusters of similar size and spread, so small subtypes are absorbed by large ones;
- k would have to be chosen, and the clusters named afterwards (the "lowest marker ratio" rule does not generalize);
- a threshold per marker can be read, checked on a histogram and changed in the dashboard, and the combinations of
  positive markers name the subtypes directly;
- a functional marker may also be expressed by tumor cells; thresholding it within the immune cells keeps those out.

Detections that are positive in too many marker channels are not cells ("Artefact"; 2.8's unfinished "crap" class):
flagged on thresholds over all cells, before the k-means, which then runs without them.

## 11. Settings and outputs

The macro saves part of its parameters in `settings.txt` (lines 85-113): `CellposeModel`,
`membraneThicknessForStraightening`, `bool_loadLabelmap` and the Cellpose environment are missing (issue I13).
The port saves the complete `settings.toml`, a `run.json` (version, counts, Cellpose version and device, timing) and
`log.txt`, plus `labelmap_final.tif` (cleaned labels) and `labelmap_classes.tif` (1 = T cell, 2 = tumor cell).

Settings file ↔ macro parameters:

| Macro | Settings |
|---|---|
| `inputImageStack`, `output` | `input.image`, `input.output` (results in `<output>/<image stem>/`) |
| `bool_limitNrCells`, `limitNrCells` | `input.max_events` |
| `ch_BF`, `ch_MOI`, `ch_tumor`, `ch_Tcell` | `channels.brightfield`, `.marker_of_interest`, `.tumor`, `.t_cell` |
| `markerLoc`, `CellOfInterest` | `channels.marker_location`, `channels.cell_of_interest` |
| `CellposeModel`, `cellDiameter` | `segmentation.model`, `segmentation.diameter` |
| `bool_loadLabelmap`, `labelmapFile` | `segmentation.method = "labelmap"`, `segmentation.labelmap_file` |
| `erosionSteps`, `minCellSize`, `maxCellSize` | `segmentation.erosion`, `.min_cell_size`, `.max_cell_size` |
| (new) | `segmentation.exclude_edge_cells` |
| `CellposeFlowThreshold` (constant) | `segmentation.flow_threshold` |
| `membraneThickness`, `membraneThicknessForStraightening` | `membrane.thickness`, `membrane.straightening_width` |
| `maxHistValue`, `binWidth` | `membrane.max_histogram_value`, `membrane.histogram_bin_width` |

---

## Issues in the macro

Kept in the port unless stated, so that the results stay comparable; to be decided one by one.

| # | Issue | Port |
|---|---|---|
| I1 | The order of the cells within an event in the ROI Manager is not deterministic (a few events on the test data). | Ordered by event, then label. |
| I2 | Cellpose runs on one montage of all events (on purpose: single events do not look like the microscopy images Cellpose is made for), as **one row**; its normalization depends on all events (so on `limitNrCells`), and a cell can cross the border between two events. | **Decided 2026-10-22:** keep a montage, but as an about-square **grid** (`segmentation.montage_columns = 0`, default); `montage_columns` = number of events gives the macro's row. Cells across a border are split (I15). Per-event segmentation was tried and dropped: slower, and further from the reference. |
| I3 | The size filter is only applied when `erosionSteps > 0`. | As the macro. |
| I4 | Tumor cluster = lowest **sum** of ratio/area (the comment says average); with unequal cluster sizes this can pick the wrong cluster. | **Decided 2026-10-22: keep the sum.** |
| I5 | `All cell combinations` gets a row only for events with cells, but the measurements are written by event index: the rows shift after an event without cells. | Fixed: one row per event. Same result when every event has cells (true for the test data). |
| I6 | ~~"Not in interface" = interface XOR membrane would also contain interface pixels deeper than the membrane band.~~ **Not an issue**: a pixel within k of the other class is within k of the edge of its own class (same k, same 4-connected steps), so the interface always lies inside the membrane band and XOR = membrane AND NOT interface. Verified on the test data: 0 pixels differ. | As the macro. |
| I7 | The erosion is from the background only: touching cells are not shrunk from each other. | As the macro. |
| I8 | The membrane analysis works on connected rings of the class mask: touching cells of the same class form one ring. | **Decided 2026-10-22: one membrane per cell.** ✅ |
| I9 | After background subtraction, pixels that are exactly 0 become NaN and can truncate a profile. | **Decided 2026-10-22: the profile length is kept explicitly; zeros are data.** ✅ |
| I10 | The control interface position is random and unseeded: not reproducible; `round(random·length − 1)` can be −1. | **Decided 2026-10-22: no fixed seed; bootstrap many random control windows per cell** to estimate the distribution of non-enriched membrane intensity: 10 windows per cell without interface, of the median interface length, pooled. ✅ |
| I11 | The average profile includes rings without an interface. | **Decided 2026-10-22: average only cells with an interface.** Verified on the reference: the macro averages all rings (its output is reproduced to 0.0005 that way); the rings without interface roughly halve the peak. ✅ |
| I12 | `Duplicate... channels=ch_MOI` on the single-channel MOI image (line 741) works only because ImageJ ignores it. | n/a |
| I13 | `settings.txt` is incomplete. | Fixed: `settings.toml`. |
| I14 | The k-means of ImageJ and scikit-learn are not bit-identical; with well-separated clusters the result is the same (100 % on the test data). | scikit-learn, seed 48, 10 initializations. |
| I15 | A cell that Cellpose draws across the border of two tiles of the montage gets one label in two events: it is measured as one cell (area of both parts) and counted in the first event. | **Fixed:** split into one cell per event, before the clean-up (`segmentation.split_across_events`). None on the test data. |
