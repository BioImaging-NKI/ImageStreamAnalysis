# ImageStream Analysis: behaviour specification for the Python port

Reference: Fiji macro **v1.8**,
`reference/fiji_v1.8/Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm`
(BioImaging-NKI/ImageStreamAnalysis, commit `cffe8fb`, 2026-10-22). Line numbers below refer to that file.
Older versions of the macro (1.5) are not used.

Validation data: a confidential test data set with the output of a v1.8 run (not in this repository; see
`tests/test_reference.py`). Numbers measured on it are kept out of this file.

Status per step: ✅ ported and validated against the reference · 🟡 ported, not validated · ⬜ not yet ported.

**Names (2026-10-31).** The port says *immune cell* where the macro v1.8 says *T cell* (as ClassifyCellTypes 2.8
does), because the immune side may be defined by several markers (`channels.immune`, summed after normalization).
Renamed outputs: `T-cell intensity` → `Immune intensity`; `Cell type` T cell → `Immune`; `T cells` → `Immune cells`;
`mean Tcell interface` / `mean Tcell NOT in interface` → `mean immune interface` / `mean immune NOT in interface`;
`T cell membrane ratio` → `Immune membrane ratio`; the table `Single Tumor Single T cell` → `Single Tumor Single
Immune cell`; `cell_of_interest` "T cell" → "Immune cell". Old settings files (`channels.t_cell`, "T cell") still load.

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

For the tumor and T-cell marker channels (and in membrane mode also the brightfield variance). The port normalizes
every immune marker of `channels.immune` this way and sums them (2.8; with one marker this is v1.8):

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
- **Reuse (2026-11-01, new):** `labelmap_cells.json` next to it holds a fingerprint (sha256 of the Cellpose input
  montage and of model, diameter, flow and cellprob thresholds, channels and the Cellpose Python) and the Cellpose
  info. With `segmentation.reuse_labels` (default on) the next run of the dataset loads the labels instead of running
  Cellpose when the fingerprint matches; clean-up, measurement and all later stages run again. This replaces the
  manual labelmap round trip for re-runs (thresholds, membrane pair, clean-up settings).

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

- `Area` = pixel count; `T-cell intensity` (port: `Immune intensity`, the summed immune markers), `Tumor intensity`
  = mean of the **normalized** marker channels.
- `Marker ratio` = T-cell (immune) / tumor intensity; `Marker ratio over Area` = ratio / area.
- k-means (ImageJ "k-means Clustering", 2 clusters, tolerance 1e-4, seed 48) on three features per cell:
  log(T-cell intensity), log(tumor intensity), area, each z-scored (mean 0, sample sd).
  Port: scikit-learn `KMeans(n_init=10, random_state=48)` (issue I14).
- Tumor cells = the cluster with the lower **sum** of `Marker ratio over Area` (issue I4).
- Class map: T cells → 1, tumor cells → 2 (`changeValues`). A third class (`class2`, green) cannot occur with 2
  clusters.
- `Class` column: in the macro the arbitrary cluster number; in the port **1 = immune cell, 0 = tumor cell**
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

`Single Tumor Single T cell` (port: `Single Tumor Single Immune cell`): the rows with 1 T cell, 1 tumor cell and a
non-NaN T-cell interface.

Plots of interface vs non-interface MOI for the cell of interest (all events, doublets): not ported (superseded by
the pair analysis, 2026-11-03). Port (2026-11-03): these columns exist once per marker of interest, with the marker's
name in brackets (`mean tumor interface (<marker>)`); see "Several markers of interest" below.

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
[−5, maxHistValue] (100 bins with the defaults) → `Histogram_*_stdDevs_above_mean.tsv/.png`. The port scales
automatically by default (see below).

## 10. Aligned profiles (lines 898-996) ✅

- Trace lengths re-measured as the first NaN of each row.
- Shift per row = 20 − (first interface pixel, or the start of an interface that wraps around x = 0) −
  interface length / 2; all four profile images rotated by it (`Array.rotate`). Rows without interface are rotated
  too (issue I11).
- `membrane_profiles_normalized.tif`; average profile over all rows up to round(mean + sd of the trace lengths) →
  `Average normalized membrane profile values - interfaces centered at 20.tsv` and its plot.

## 8-10 in the port (`membrane.py`)

Optional (`membrane.analysis`, default off: not every ImageStream analysis needs it).

- **One membrane per cell (I8)**: the outline of the cell (marching squares at 0.5 on its own mask, i.e.
  along the pixel edges; smoothed along the line with a Gaussian, σ = 2 samples, in place of Fiji's spline fit;
  resampled at 1 px), straightened with `membrane.straightening_width` samples across it (bilinear), maximum
  across → profile. Background: the MOI median of the stack. The profile lengths are kept (I9).
- **Any pair of cell types (2026-11-01).** Every cell but the artefacts gets its profile and a *neighbour track*:
  for each position, the label of the nearest other cell within `membrane.thickness` (4-connected steps, the reach
  of the interfaces) of the cell's own pixels there, straightened along the same line (nearest sample, maximum
  across), 0 = none. `membrane.cells` and `membrane.partners` choose the pair (a cell type, an immune subtype
  "Immune A+ B+", or a marker "A+"; empty: the cell type that carries the marker, from the data (2026-11-03),
  against the other type). Interface = the
  positions whose neighbour is a partner. **Non-interface = the membrane that touches no cell** (v1.8: all of the
  membrane outside the interface, so also contacts with cells of the same type); the control windows are taken on
  that free membrane too, read as one closed line. Reason: a contact with another cell that carries the marker
  (e.g. two tumor cells) is brighter (two membranes side by side) and would raise the baseline, and controls
  that fall on it looked enriched. On the test data this raises the median `stdDevs above mean` by about 0.2 SD;
  the control quartiles stay within 0.15 SD of Fiji's. (With the v1.8 baseline the neighbour tracks reproduce the
  previous port to a few hundredths.)
- **The band centred on the membrane (the user, 2026-11-04; `membrane.band_centre`, default "membrane").** The
  outline does not lie on the membrane: at a free membrane Cellpose's mask reaches beyond it, at a contact the
  outline is the boundary between two masks. A band centred on the outline (v1.8; "outline") reaches less of the
  free membrane than of an interface, and a narrow band fakes enrichment (the user test of 2026-11-01 led to the
  4.5 µm default). Now `membrane.locate` first samples every marker of interest across the outline (±3 µm, 1/4 px
  steps, outward normal from the own mask), averages it per cell along the free membrane and along the contacts with
  each cell type (positions as the neighbour track, ≥ 5 per cell), then over the cells of each type (≥ 20). The
  peak of the average (parabolic refinement; none at the window edge) is the membrane's place per (cell type,
  place). The band's centre line is the outline moved along its normal by that offset (per position by its place,
  smoothed over 2·reach+1 samples); the neighbour track stays on the outline. `Profiles.lines` is that centre line
  (the dashboard's band). Choices:
  - **Per cell type and place, not per cell** (the user, 2026-11-04): robust, and the marker does not choose its own
    band cell by cell. On the real data the between-cell spread was larger than the within-cell spread (local
    notes), so a per-cell offset remains an option to examine.
  - **Only markers that the cell type carries place its band** (`carrier` from whole-cell intensities, before the
    membranes). On real data, a marker on the other type showed no membrane peak in a cell type, or a peak inside
    the cell (intracellular signal); the clearest peak (prominence: peak above the window ends over the range, ≥ 0.3)
    among the carried markers places the band, and holds for every marker. No carried marker, too few cells, no
    peak: on the outline, logged.
  - **The width stays** (about 4.5 µm): centring removes the systematic free-vs-interface difference, the width
    takes up the cell-to-cell variation.
  Findings on real data (numbers in the local notes): the outline lies outside the membrane by a similar distance at free
  membranes of both types; for one type (the immune cells) also at contacts, so the outline-centred band was fair for
  it; for the other (the tumor cells) the membrane at contacts with immune cells lies on the partner's side of the
  mask boundary, so the outline-centred band favoured the interface. Consequence for `inside_only`: it cuts away part
  of the measured cell's own membrane at such contacts (a further reason it is only a sensitivity test).
  Checks: the log per (type, place) with the spread of single cells; `Membrane location.tsv` and `Membrane location
  - profiles.tsv`; the dashboard's *Across the membrane* (membrane card) with the band of each place; with "outline"
  a WARNING when the band's outermost sample is less than a pixel beyond a membrane. Synthetic test
  (`test_the_band_follows_the_membrane_not_the_outline`): an evenly stained membrane 2 px inside a free outline and
  on the contact line; 3 px band on the outline: interface/free ratio 1.3, centred: 0.93 at 3 and 9 px. The Fiji
  comparison (`test_reference.py`) keeps "outline".
- **Where the marker of interest is expressed (the user, 2026-11-01).** If the partners carry the marker too,
  their membrane lies against the measured one at the interface (half the straightened band lies in the partner)
  and looks like enrichment even when the staining is even. `expression_check` compares the median free-membrane
  intensity of the partners with that of the measured cells (without profiles: the whole-cell means): ≥ 1 is a
  WARNING (expressed in the partners? check `membrane.cells` / `membrane.partners`), ≥ 0.25 a CAUTION.
  In the log, `run.json` ("marker of interest") and the overview (`MOI partners / measured`).
  `membrane.inside_only` (default off, as the macro) takes the maximum across the band only over the samples inside
  the measured cell. On the test data the check says CAUTION, and inside_only lowers the median `stdDevs above
  mean` and the median interface/non-interface ratio clearly, while enrichment remains (numbers in the local notes).
  **Caveat (the user, 2026-11-01):** touching membranes are nanometres apart, far below the resolution of the
  ImageStream, so the blur mixes both membranes in every pixel near the interface; inside_only changes which pixels
  are sampled, not what they contain, and cannot separate the two. The expression check is the safeguard; the
  option is kept as a sensitivity test.
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
- **Automatic histogram scale (2026-11-02, the user).** `membrane.max_histogram_value` and
  `membrane.histogram_bin_width` default to empty = automatic (the macro: 20 and 0.25, 100 bins).
  `membrane.histogram_scale`: maximum = the 99th percentile (linear interpolation) of the finite interface values or
  of the control windows, whichever is higher, at least 5 (20 without values), rounded up to a whole number of bins
  from −5; bin width = the smallest of 0.1, 0.2, 0.25, 0.5, 1, 2.5, 5, 10, 25, ... giving at most 60 bins (all divide
  5, so a bin edge lies at 0). Explicit values are used as before (round((max + 5) / width) bins). Reason: with a
  fixed 20 SD most of the plot was empty or, with strong enrichment, the tail was cut; a percentile keeps one
  extreme cell from squeezing the rest. The lower end stays at −5, as in the macro (the control windows centre on 0, and
  strongly negative interfaces are rare). Effect: only the binning of the two histogram tables changes; no cell
  value changes (not rerun on the test data). The dashboard repeats the rule in JS (`histScale`), with −/+ buttons for the maximum and the bin width and
  *auto*, rebinning live from every cell's value and the control windows; the scale is kept in the remembered view.
  Its control windows use a fixed seed, so its automatic scale can differ by a bin from the run's tables.
- `membrane_profiles_normalized.tif`: one row per row of the membrane table, NaN beyond the length. Plots ⬜ (dashboard).

Validation on the test data (2026-10-22), on distributions: the number of cells with an interface and the medians
of the membrane length, the interface and non-interface means and their ratio are within 1-10 % of Fiji's; the
quartiles of `stdDevs above mean` within 0.5 SD; the pooled control windows have the quartiles of Fiji's single
controls within 0.25 SD; the average profile of the cells with an interface correlates > 0.99 with Fiji's profiles
of the rings with an interface. (These are the checks of `tests/test_reference.py`.)

## Projects and the import page (new) ✅

- `project.py`: a project = several IDEAS zips (of one or more experiments) with shared analysis settings, saved as
  `project.toml` in the project folder. The roles (brightfield, marker of interest, tumor, immune and functional
  markers) form a **panel by marker name** in the project's `[channels]` (2026-10-31); per zip: path, pixel size, the
  channels to load and their names, and per dataset an experiment label and whether to include it. Each zip's names
  give its channel numbers, so exports with differently numbered channels share one panel. Projects written before
  had roles by number per zip (`t_cell` among them); on loading they become the panel by name when all zips agree,
  and otherwise stay as per-zip overrides. On the page a role is a tick box per named channel, and holds for every zip
  with a channel of that name.
  `imagestream-cli run-project` analyses every included dataset into `<project>/<experiment>/<dataset>/`; a failing
  dataset is logged (and recorded in its `summary.json`) without stopping the others. `overview.tsv` holds one row
  per dataset (counts, doublets, median membrane ratio of the doublets, and the membrane results if computed).
- `gui.py` + `templates/gui.html` (`imagestream-cli gui`, or `imagestream`): a local page to add zips, name the
  channels with example images of each, give them a role, set the pixel size and the analysis settings, run with the
  log streamed, and see the overview. Server: stdlib `ThreadingHTTPServer` on 127.0.0.1 with a per-session token;
  stops after 15 min without requests; a second start reuses the running one.
- **Run state and stale results (2026-11-04).** The user found the overview showing the previous results during a
  run, with no clear sign of a run going on, and old results lingering in a reused output folder (the overview reads
  every `*/*/summary.json` of the folder).
  - The page shows the run's state under *Run all datasets* and next to the title of step 5: *Running…* (spinner,
    progress bar, dataset k of n, stage with pools per experiment, time), then *Done*, *Done, N of M failed*,
    *Stopped* or *Failed (exit code N)*. The start buttons are off while it runs; *Run all datasets* reads *Running…*.
  - `gui.Job` records the result folders a run writes (`run_folders`, from the project and `--dataset`) and its start
    time. A row of the overview is *running* (in the run, `summary.json` older than the start), *updated* (newer),
    *previous run* (the run ended before it: stopped or failed) or *waiting* (no folder yet). The overview is refreshed
    on every progress line and at the end; the partner preference and the comparison are greyed out during a run
    (they are written at the end).
  - **Results not part of the project** (`project.in_project`: the folder is not `<experiment>/<dataset>` of an included
    dataset: removed, unticked, or experiment renamed) are left out of `overview.tsv`, `interactions.tsv`,
    `subset preference.tsv`, `membrane pairs.tsv` and the comparison, and hidden on the page behind *Show them*.
    Before, they were pooled into the partner preference of their experiment. Their files are never deleted.

## Fixes from the user test (2026-11-01)

Four agents used the tool from the manual alone, on synthetic exports with a known truth (scenarios outside the
repository). Fixed from their reports:
- **No marker of interest needed** for cell types and contacts (`channels.marker_of_interest = []`): the event table
  leaves out the marker columns and gets `Tumor-immune interface` (doublets use it), no marker check, the dashboard
  hides the marker cards; the membrane analysis then refuses to run.
- **Marker check**: partners / measured 0.8–1.25 = "on both at a similar level" (membrane doubling; compare with
  same-type contacts), ≥ 1.25 = brighter on the partners; measured cells below 25 % of the brightest cell type =
  NOT INFORMATIVE (grey badge), unless the partners are brighter.
- Overview: `membrane pair` column; subtype columns 0 instead of empty, in one order; no `inf`; hover texts.
- Dashboard: colour slots for the immune types present (not for all 2^n combinations); no `getImageData` on a canvas
  of width 0; role names as on the page; the subtype card says the threshold is an intensity.
- Import page: an unsaved setup is kept as a draft in the browser until saved.
- Log: "subtype markers" (not "functional"); with one subtype marker the artefact rule is off and the log says so.
- Manual: the control windows are no control for the interface (compare pairs; geometry alone gave ~1.3× / 1–1.5 SD
  in the test), thresholds can be off both ways, sizes in pixels, overview columns, comparing conditions.

Decided by the user and built (2026-11-01, same day):
- **Threshold "Two populations" (default)**: noise limit of the negative cells (KDE mode + 2.5 robust SD from the
  lower side) and a 2-component Gaussian mixture; threshold at the crossing of the
  weighted densities when the upper component is credible (≥ 0.5 %, Ashman's D ≥ 2, mean above the noise limit), else
  the noise limit. Simulation (500 samples): errors 4790 vs best 2417, MaxEntropy strict 14 881, MaxEntropy 21 965.
  The user: positive peaks are often small, which is why a stricter threshold worked better.
- **Straightening width default about 4.5 µm** (`membrane.straightening_width` empty): 5 px at 40× gave 1.4–1.6 SD of
  false enrichment at contacts without biology (outline outside the membrane at free edges); 9 px 0.35 SD, real
  enrichment kept. A reference channel was tested (the cells' own cytoplasmic tracker) and rejected: 4.8 SD false.
- **Pair table** in the dashboard: median enrichment for every Membrane of × Interface with, with the marker check.
- **Subset preference** (`interactions.subset_preference`, within events with ≥ 2 immune cells, bootstrap 95 %
  interval, 1000 resamples, seed 48; page 400): `Subset preference.tsv`, project `subset preference.tsv`, dashboard
  tab "among immune cells". The dataset-wide variant was tested and is too insensitive for conjugate exports.

## Fixes from the second user test (2026-11-01)

Round 2: three of four users reached the right conclusion; counts matched the truth almost cell for cell. Fixed:
- **Artefact rule** on the subtype thresholds of the immune cells of a first k-means (rows "artefact rule" in
  Marker thresholds.tsv), then the final k-means without the artefacts. Over all detections the thresholds lay at the
  noise: 5 of 47 "artefacts" were cells with a neighbour's light, and the user then switched the rule off, which made
  the 42 real artefacts a fake triple-positive subset. Now exactly the 42.
- Subset preference: no interval below 5 expected contacts (zero-width intervals coloured 1-contact squares).
- `Membrane pairs.tsv` / project `membrane pairs.tsv`: every pair (membrane.pair_table), so no re-run per pair.
- Overview `marker check` in words (`expression_check(...)["level"]`); the macro's doublet ratio and the per-event
  card carry a caveat.
- `membrane.thickness` default about 1.5 µm (2/3/5 px at 20/40/60x); the reference comparison keeps its 2 px.
- Draft of the import page also saved on pagehide; favicon; no tracebacks for aborted connections; one subtype
  marker: no duplicate marker selection; artefacts counted on the cell card; warning above six immune subtypes.

## Comparing conditions (2026-11-01) ✅

`comparison.py`; per dataset `condition` and `replicate` in project.toml, `[comparison] reference`. Units =
replicate × condition (datasets pooled). Subset preference: y = ln((O+.5)/(E+.5)), tested when E >= 5 in every unit;
membrane enrichment: cell-weighted median SD, pairs of the cell type that carries the marker with the other type
(per marker of interest, from the column `marker on`; majority of the datasets, with a note, 2026-11-03), >= 10
interfaces.
Paired t over replicates with both conditions (>= 2), else Welch (>= 2 per side), else no test; BH q per measure and
condition. Written to `condition comparison.tsv` (run end, Save, `imagestream-cli overview`); the overview computes it
live from the page's labels. Decided **no test per donor** (two samples of one donor against each other with the
event-sampling error): calibrated only without any sample-to-sample variation, 36-75 % false positives with an SD of
0.03-0.1 on y. Calibration: `tools/stats/calibrate_comparison.py` (results in the manual): paired 4.3-5.9 %, Welch
2.4-4.7 %, BH any-false 4.9-5.8 %; power for the user-test effect (true 1.07 -> 1.12) 11 % with 3 donors x 300
events, 74 % with 4 x 3000. The jackknife over events matched the true SD of y (the linearized SE was 8 % low at 150
events), but is not used, since there is no per-donor test.

## The dashboard (new; replaces the macro's plots and overlays) ✅

`dashboard.py` + `templates/dashboard.html`: every analysed dataset gets `dashboard.html` (results embedded as JSON)
and `dashboard_files/` (event images: per channel, 8-bit, display range background (median) to the 99.8th
percentile, brightfield 0.2-99.8th percentile; sheets of 10 x 10 events, WebP). It opens from disk (file://) or
from the import page (served from the project folder, token-protected). No external resources.

- Gallery of the events (only the visible tiles are drawn): channels composed in the page (colour per channel,
  on/off, brightness; brightfield starts at half brightness), cell outlines in the colour of their type,
  interfaces, and optionally *contacts*: a line between the centres of touching cells.
- **Cell types (2026-11-01)**: tumor, immune (or one type per immune subtype), artefact. Colours follow the entity in
  a fixed order (dataviz palette): tumor slot 2, the immune types slots 1, 4, 5, 6, 7, 8 in the order of their
  combination of positive markers (bit order of `channels.functional`, "negative" first), more than six: grey
  "other"; artefacts muted grey; the interface slot 3. Legends name every type; in the class scatter (more than
  three colours fails the all-pairs colour-blind check) a click on a legend entry brings that type forward.
- **Subtype markers (live)**: per functional marker a histogram of log(v − min + 1) of the immune cells with the
  run's threshold (dashed) and a slider; moving it re-labels the subtypes everywhere (gallery, plots, contacts,
  membrane pairs, tables). The tumor / immune split and the artefacts stay as computed. "Copy the thresholds as
  settings" gives `[classification] thresholds = {...}` and `[membrane] cells / partners` to paste. **Apply to the
  project** (2026-11-01; only when served by the import page, which has the token in the URL `/files/<token>/...`):
  the moved thresholds (null = automatic again) and a changed membrane pair go into `project.toml`
  (`gui.apply_from_dashboard`, project-wide), the project runs again (Cellpose reused), the page reloads when done.
  The import page notices a run started this way and reads the project settings again from disk.
- **View** (2026-11-01): channels, overlays, filter, sort, size, pair and tabs in localStorage per dashboard path
  (the server token stripped). **Download .tsv** of the current table (all rows, filtered and sorted, 3 decimals).
- **Interactions**: a contact matrix (cell type × the type it touches), coloured by observed / expected (diverging,
  ¼× to 4×, log scale; expected exactly as `interactions.contact_table`, so page and TSV agree) or by the fraction
  touching (sequential); a click on a square filters the gallery to the events with that contact (a chip in the
  filter row clears it) and draws those contacts. Partner bars: per cell type, its contacts by neighbour type
  (share or number), and the share of its cells that touch any cell. A Contacts tab in the table.
- Cell classes: area against marker ratio, log-log (the macro's `cell Classes.png`).
- (Removed 2026-11-03.) The per-event plot of the MOI at the interface against the rest of the membrane (the
  macro's `Channel N in ... membrane plot`); its values stay in the Events table and the sort.
- **Membrane analysis (any pair, live)**: "membrane of" and "interface with" (a cell type, an immune subtype, or a
  marker); the page computes the enrichment from every cell's profile and neighbour track as `membrane.py` does
  (the control windows with a seeded generator, so they differ from the TSV's random ones; on the test data the
  medians agree to 0.01 SD), with the histograms, the average aligned profile and the heatmap, and the expression
  check as a badge (icon and words: ok / partners carry part of the marker / marker brighter on the partners).
- Cell card (with the subtype markers, the touching cells and the membrane), cell, event and contact tables
  (sortable). Hover and selection are linked in all views; one filter row (events: all / both cell types / one
  tumor + one immune cell / with an interface / with touching cells; a contact; sort) scopes everything. `#cell=N`
  or `#event=N` in the address selects on opening.
- **Cell card: channel montage and membrane band (the user, 2026-11-04).** Under the measurements, every loaded
  channel of the event in grey (its display range, also when switched off in *Channels*), with the selected cell's
  outline. *membrane band* (space bar; remembered in the view) lays the band that was straightened over the image and
  the montage at 40 % opacity: `Profiles.lines` (the smoothed outline of `membrane.outline`, one point per profile
  sample, in event pixels) are written to the page as `membranes.line` (int16, 1/8 px) with `band_width` (the
  straightening width in px), and the page draws `band_width` / 2 on each side of the normal as `straighten` samples
  it (with `inside_only` only the inner half). Colour per position, from the neighbour track and the pair shown:
  interface with a partner (aqua), free membrane, the reference (amber), another contact, left out (grey). A selected
  cell of the partner type shows its contacts with the measured cells as interface; a selected event shows the measured
  cells of the pair only. Size: 4 bytes per membrane sample (the profiles take 8 per marker).
- Colours: the validated reference palette (dataviz skill), in light and dark mode; aqua and yellow are below 3:1
  on the light surface, so they are always named in a legend and the values are in the tables.
- Project overview (import page): the contact matrix of each experiment against chance, side by side (from
  `interactions.tsv`).
- Cost on the full test dataset: about 6 s; outlines are traced on each cell's bounding box.

## Several markers of interest; the measured cells from the data (the user, 2026-11-03) ✅

Decisions by the user, 2026-11-03, and their reasons:

- **Enrichment for any marker.** `channels.marker_of_interest` is a list (old files: a number or a name becomes
  `[x]`, `""` becomes `[]`). Every marker gets the whole membrane analysis (profiles, `analyze_membranes`,
  `pair_table`, `expression_check`). The outlines and neighbour tracks are computed once
  (`membrane.membrane_profiles_of`; the tracks are shared objects); only the profiles differ per marker.
- **No compatibility with the macro's names needed** ("as long as it's scientifically sound"). There is no primary
  marker: every marker is treated alike, also when there is only one. Membrane files:
  `<macro name> - <marker>.<ext>` (`Membrane intensity analysis - <marker>.tsv`, the average profile, both
  histograms, the control windows, `membrane_profiles_normalized - <marker>.tif`; the name made safe for the file
  system; a name that two markers share gets ` (ch<N>)`). `Membrane pairs.tsv` stays one file with the columns
  `marker` and `marker on` in front. `run.json` lists every marker ("markers of interest": carrier, medians, pair,
  marker check, membranes); the former single "marker of interest" entry is gone. Overview: the columns of each
  marker end in ` (<marker>)` (`MOI on`, `MOI partners / measured`, `marker check`, `median <type> membrane ratio,
  doublets`, `membrane pair`, `membranes`, `with interface`, `median stdDevs above mean`, `median control (windows)`).
- **Per-event columns kept, per marker.** `All cell combinations.tsv` and the doublet table have the six MOI columns
  once per marker, as `<macro column> (<marker>)`. Reason: they cost nothing, keep the Fiji comparison of the event
  tables possible (the reference test maps the names back), and stay useful as a quick look; but they are the
  macro's measure (all partners together, no control), so nothing else depends on them.
- **The measured cells follow from the data.** `channels.cell_of_interest` is removed (an old file loads; the value
  is dropped with an info line). Empty `membrane.cells` means, per marker, the cell type (Tumor or Immune) with the
  higher median of the measure of `expression_check` (free membrane, background subtracted; without profiles the
  whole-cell mean): `membrane.carrier`. Empty `membrane.partners`: the other type. Explicit cells / partners apply to
  every marker. Recorded in the log ("<marker>: carried by the tumor cells (medians ...): membranes of Tumor with
  Immune"), in `run.json`, in `Membrane pairs.tsv` (`marker on`) and in the overview (`MOI on (<marker>)`). Reason:
  a wrong setting measured the partners' membrane; the data say where the marker is. On the test data the carrier
  is the macro's CellOfInterest (checked in `test_reference.py` from the macro's `settings.txt`), so the
  reference comparison is unchanged. Consequence: with the default pair the marker check cannot say "brighter on the
  partners".
- **Comparison of conditions:** one measure per marker, `membrane enrichment (<marker>)`, on the side of the type
  that carries it (column `marker on` of each dataset; when the datasets disagree, the majority, with a note on its
  rows). Results without the column (before 2026-11-03) are not compared, with a note; the statistics are unchanged
  (no recalibration needed: the test itself did not change, q stays per measure).
- **The dashboard's per-event panel ("Marker of interest at the interface") is removed:** superseded by the pair
  analysis (all partners together and no control made it misleading). Its values stay in the Events table (per
  marker) and in the sort "membrane ratio per event" (of the marker shown). The membrane card has a *Marker* menu
  (two or more markers); the tracks are embedded once, a float32 profile per marker; each marker keeps its own pair
  on the page (default: that of the run); the chosen marker and the pairs are remembered in the view. Apply /
  Copy as settings write one pair, which applies to every marker (said in a tooltip and the exported text).
  Superseded the same day by a pair per marker (below).
- **Import page:** the MOI role takes several channels; *Marker of interest expressed in* is gone; *Membrane of* /
  *Interface with* are menus (automatic first, stored as `[]`; Tumor, Immune, every combination of the ticked subtype
  markers as `Immune A+ B+` / `Immune negative`, and with two or more markers `A+`: the strings of
  `membrane.select` and of the dashboard's `selectorOptions`; the dashboard lists only the combinations present).
  A stored value not offered stays as its own option "(from the settings file)".
- **The pair is validated** (a bug the user hit: `membrane.cells = ['Immune A']`, typed by hand, matched no
  cell silently, and saving an empty profile stack crashed). `membrane.selector_problem`: after "Immune " come
  subtype markers of the panel as "<m>+" in the panel's order, or "negative"; otherwise an error with the likely
  meaning ("did you mean 'A+' ... or 'Immune A+' ...?"), in `pipeline.start` (before segmentation) and in
  `Project.problems()` (the import page lists it too). A valid subset without cells in a dataset only warns; that
  marker then writes no membrane files (and removes those of an earlier run), and the dataset finishes normally.
- Old project folders: the overview reads them (their `Membrane pairs.tsv` gets an empty `marker`, with a log
  line); the comparison skips their membranes with a note.

## A membrane pair per marker of interest (the user, 2026-11-03) ✅

Reason: a pair set for all markers did not fit a panel with one marker on the tumor cells and one on an immune subset
(the first measured on the tumor cells against the subset, the second on the subset against the tumor cells).

- `membrane.pairs`: `{ <marker name> = { cells = [...], partners = [...] } }`, written by `Settings.to_toml` as
  sub-tables `[membrane.pairs.'<marker>']` (so that a table pasted from the dashboard fits a file; an inline table
  could not be extended). `MembraneSettings.migrate` checks it (only `cells` / `partners`; a string is a list of one).
- Resolution per marker and per side (`MembraneSettings.pair_settings`, `membrane.pair_of`): the entry in `pairs`,
  else the global `membrane.cells` / `partners` (kept: compatibility, and one pair for all markers in a settings
  file), else automatic (the carrier, then the other type of the carrier).
- Checks (`selector_problems`, in `pipeline.start` and `Project.problems()`): the selectors of every entry, and an
  entry for a name that is not a marker of interest.
- `pair_table` adds the run pair of each marker when it is not one of the table's rows × columns (a marker such as
  `A+` as partner, or several selections), so `Membrane pairs.tsv` and the comparison of conditions have it. The
  comparison still takes the side from `marker on` and the partners from the other type; an explicit pair of the
  carrier against the other type (or a subset of it) is therefore compared, a pair within one type is not.
- Import page: one row per marker of interest ("<marker>: membrane of [menu] interface with [menu]"), no global
  row. The rows follow the MOI ticks and renames (a renamed channel keeps its pair; an unticked marker loses its
  entry). An older project's global pair is turned into an entry per marker when the page loads the project
  (`Project.pairs_per_marker` in `gui.project_to_json`) and saved that way; the file itself changes only on saving.
- Dashboard: *Apply to the project* sends `{marker, cells, partners}` and `apply_from_dashboard` writes
  `membrane.pairs[<marker>]` only; *Copy as settings* gives `[membrane.pairs."<marker>"]`.

## The excess over both membranes (the user, 2026-11-03) ✅

Question: correct the enrichment for a marker that the partner carries too. Assumption: at an interface the band
holds both membranes (closer than the resolution) and background-subtracted intensities add. Per cell with an
interface: `excess = (mean(p[interface]) - own free mean - mean over the interface positions of
partner_free[track]) / own free sd` = `stdDevs above mean - mean(partner_free) / own free sd`, where `partner_free`
is the partner's mean free membrane (`membrane.free_means`, the marker check's per-cell value).

- Built side by side with today's score: the column `excess stdDevs above both membranes` (membrane table), the
  aligned excess profile (`MembraneResult.excess_profiles`; the average as `excess profile`), and
  `median excess above both membranes` per pair in `Membrane pairs.tsv`. The dashboard has a *Measure* menu ("SD
  above own membrane" / "SD above both membranes") that the histograms, average, heatmap, cell card, *All pairs*
  and the Cells table follow; remembered in the view. JS parity with Python on the two-marker synthetic run: within
  the 3 decimals of the TSV (scores, average excess profile, pair medians).
- Choices: a partner position whose partner has no free membrane is left out of the partner mean (the aligned excess
  profile subtracts the cell's partner mean there); no partner with a free membrane: the excess is NaN. The control
  windows have no partner, so their excess is the control: no separate column.
- Simulation (`tools/stats/simulate_excess.py`, `test_excess_over_both_membranes_by_simulation`): thin membranes
  drawn on a 4x finer grid (so that touching membranes coincide), blur sigma 1.2 px, width 9 px, reach 3 px.
  Medians (30 events each): (a) marker on both, uniform: score 7.6, excess -1.0; (b) marker on the measured cell
  only, enriched 3x at the contact: score 14.9, excess 14.0; (c) on both, enriched 3x on the measured side: score
  22.4, excess 13.8. With membranes drawn 2 px wide each (side by side, so 2 px apart), (a) gave excess -6.0: the
  maximum across the band is additive only when the two membranes coincide within the blur.
- **On real data the assumption failed** (numbers in the local notes): contacts between cells of the same type,
  where no enrichment is expected, gave a clearly negative median excess, for both markers of interest. At the
  ImageStream's resolution and with the maximum projection, the partner's free membrane over-corrects. The
  measure is therefore offered as a conservative check (it can hide an enrichment, not invent one), with the
  advice in the manual to check same-type contacts first.
- **Comparison of conditions:** stays on `stdDevs above mean`. The excess should not become a second measure as
  long as same-type contacts do not give about 0 on real data: a second test of a measure that is biased would add
  multiple testing without a clear meaning. Revisit if a better correction (for example a scale factor fitted on
  same-type contacts, or the sum instead of the maximum across the band) makes same-type contacts about 0.
- **The sum across the band (the user, 2026-11-03).** `Profiles.sums` holds every profile as the sum over the
  straightened band next to the maximum (`Profiles.summed()` swaps them; with `inside_only`, the samples outside the
  cell count 0). The dashboard has a switch *across the band: maximum / sum* that the whole membrane card follows
  (both measures, *All pairs*, the marker check), remembered in the view; the tables and `Membrane pairs.tsv` keep
  the maximum. Simulation: (a) excess +0.9 with the sum (-1.0 with the maximum), (b) and (c) about the same; the
  test checks both projections. On the real dataset (numbers in the local notes) the sum brings the same-type
  contacts clearly closer to 0 for both markers, but three of the four stay negative: better, not yet sound. The
  excess therefore stays experimental and out of the comparison of conditions.
- **A scale factor: examined, not adopted (the user, 2026-11-03; may be revisited).** Idea: subtract *k* times the
  partner's free membrane, with *k* fitted so that same-type contacts give 0. Diagnostic on the real dataset (script
  and plot kept locally): per cell with a same-type interface, y = interface mean - own free mean against x = the
  partner's free mean over the interface positions.
  - Maximum across the band: not linear (flat for dim partners, then a slope well below 1). No factor can fix it.
  - Sum across the band: about linear and through 0, for both markers, so the partner's level is taken into account
    as the model assumes. But the slope is below 1 and depends on the cell type of the partner, not on the marker:
    clearly lower at immune-immune than at tumor-tumor contacts, alike for both markers. That points to geometry (the
    share of the partner's membrane that falls in the band; small, curved cells give less). For a marker on both
    types, the slopes of the mixed contacts were close to the same-type slope of the partner type.
  - Per cell the scatter is large (correlations about 0.6-0.8): the correction holds on average, not per cell.
    The own free membrane as a second regressor has a small negative weight (the own membrane counts somewhat less
    than fully at an interface).
  - Not adopted because: fitting on same-type contacts assumes they are not enriched (doubtful for a homophilic
    adhesion molecule, where it would remove real enrichment); the mixed pairs, the ones of interest, have no neutral
    reference to check a factor per partner type; and two fitted parameters per dataset with a hidden assumption
    are hard to defend. If revisited: a factor per partner type, fitted on same-type contacts with a marker that is
    not expected to be enriched there, and checked on a dataset it was not fitted on.
  - Both the interface and the free membrane are means *along* the membrane (and the SD is over positions), so the
    sum *across* the band (a fixed number of samples) does not make the measure depend on the cell's size.

## Zero intensities in the clustering (2026-11-03)

One IDEAS export had a block of zero pixels in one channel of one object; a cell inside it had a mean tumor intensity
of exactly 0, so log(0) stopped the k-means ("Non-finite features"). `kmeans_classes` now gives a cell with an
intensity of 0 or less half the lowest positive intensity of that marker, and logs its events; `assign_cell_types`
skips the infinite marker ratio of such a cell in its cluster sums. In that dataset the cell then clusters as an
immune cell (its immune intensity is high); nothing changes for data without such cells. The macro would have
clustered a -Infinity feature (ImageJ k-means behaviour undefined).

## Merging ClassifyCellTypes 2.8 (decided 2026-10-31; classification and interactions built)

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

**As built (`classification.py`, 2026-10-31) ✅ (synthetic tests; thresholds equal Fiji)**

1. *Thresholds* (2.8, lines 929-1000): per functional marker, on log(v − min + 1) of the background-subtracted
   cell means (`Mean ch<N>`), with ImageJ's auto-threshold on the 8-bit histogram, as for the Otsu normalization:
   `classification.threshold_method` "MaxEntropy" (2.8's choice), "Otsu", "Yen" or "RenyiEntropy" (all four ported
   from ImageJ's AutoThresholder and equal to Fiji), or "MaxEntropy strict" (new, 2026-10-31: the MaxEntropy
   threshold, raised if needed to the median + 2 robust SDs (1.4826 × MAD) of the values below it, on the log
   scale); back-transformed, exp(t) − 1 + min. On synthetic log-normal mixtures Yen and RenyiEntropy were not
   stricter than MaxEntropy, and Otsu was more liberal when positives are rare; "MaxEntropy strict" roughly halved
   the false positives of well-separated markers for a few more misses.
   Positive: v > threshold. MaxEntropy is transcribed from ImageJ's AutoThresholder (Kapur, Sahoo and Wong);
   verified against Fiji: equal to 1e-11 on 12 log-normal mixtures, and a fixed case is in the tests. A threshold
   can be set by hand per marker name (`classification.thresholds`). MaxEntropy cuts in the upper tail of the
   negative population, so it calls a few negative cells positive.
2. *Artefacts*: a detection positive in at least `classification.artefact_min_positive` functional markers, with
   the thresholds over **all** detections. Default (the user, 2026-10-31): positive in all functional markers, when
   there are two or more; 0 = no artefacts. Cell type `Artefact`, Class 2, class map 3;
   not counted as tumor or immune cells, not part of the interfaces. Only the functional markers count. (Counting
   the tumor and immune markers too, so that a detection bright in both would be an artefact, was tried and dropped:
   it needs a positive/negative split of those markers, and on the test data an automatic threshold over all
   detections calls far more cells positive than the k-means finds. Tumor against immune stays the k-means alone.) 2.8 had: positive in all functional
   markers (`removeBadDetections`, switched off, and it read the positives before they were computed).
3. *Lineage*: the k-means of v1.8 on the other detections (features and the tumor-cluster rule unchanged). 2.8
   clustered on the sum of the phenotypic markers, the tumor marker, optionally the functional markers (added
   unnormalized: a mistake, line 650) and area, and picked the tumor cluster by the sum of the ratio; not taken over.
4. *Subtypes*: the thresholds of step 1 again, now computed on the **immune cells** only; columns `<marker>+` (0/1)
   and `Subtype` of the immune cells: the positive markers ("A+ B+"), or "negative". `Positive markers` holds the
   count of step 2.
5. *Pools* (`classification.pool`, projects only): "experiment" (default) classifies the cells of all datasets with
   the same experiment label together: one set of thresholds and one clustering, so subtypes compare between
   datasets. The markers are matched by role and position, so the datasets may number their channels differently.
   `run_project` then measures every dataset first, classifies each experiment, and finishes every dataset (the
   images are read again; the log continues). "dataset": each on its own, as `run`.
6. Outputs: `Marker thresholds.tsv` (marker, role, cells it was computed on, method, threshold, positive, of);
   `run.json` and the overview count artefacts and subtypes and say what a dataset was classified with. The
   dashboard draws artefacts in grey and shows the subtype (the colours per subtype come with the dashboard work).

**Interactions (`interactions.py`, 2026-11-01) ✅ (synthetic tests; consistent with the interfaces on the test data)**

- 2.8 records which cells touch with MorphoLibJ's Region Adjacency Graph on the eroded montage labels (lines
  479-484, 609-632): `Interaction with (IDs)` (1-based row numbers), the subtypes of up to four partners, the number
  of interactions per image, and a `multiplets` table; the notebook `Summarize_Interactions.ipynb` then
  cross-tabulates the partner types.
- The port: two cells of the same event interact when one comes within `membrane.thickness` pixels (4-connected
  steps) of the other, the reach of the interfaces, so a cell with an interface always has its partner in the list
  (verified on the test data: every event with an interface has an interaction). Artefacts take part as detections.
  Cells are named by their label (`Label`, unique in the dataset), not by row number.
- Cell table: `Neighbours`, `Touching (labels)`, `Neighbour types` (each partner's type: Tumor, Artefact, Immune or
  "Immune <subtype>"). Event table: `Artefacts`, `Interactions` (pairs), at the end, after the macro's columns.
- `Interactions.tsv` (replaces the notebook): per cell type and neighbour type, the cells of the type, the contacts
  (each pair counts from both sides), the cells touching at least one such neighbour, and that fraction (all
  events); and the **partner preference**: `contacts (3+ cells)` against `expected (3+ cells)` and their ratio.
  `Multiplets.tsv`: the number of events with n cells and with n interactions.
- **Partner preference (decided with the user, 2026-11-01).** Expected = the mean over all ways of shuffling the cell
  types among the cells of each event (exact; same in the page). Only events of three or more cells: in an event of
  two touching cells the pair can only be what it is, observed = expected by construction, and such events pulled
  every ratio towards 1 (on the test data about half of all contacts). Events without contacts, and the events left
  out in IDEAS, take no part: it is a preference among the exported multiplets, not a binding frequency (that needs
  an export with the singlets). Called so in the page, the overview and the manual.
- Project: `interactions.tsv` with the tables of all datasets and their sum per experiment (dataset "(all)"; the
  ratio from the summed counts); the overview gets `interactions` and `tumor–immune contacts`.

## 11. Settings and outputs

The macro saves part of its parameters in `settings.txt` (lines 85-113): `CellposeModel`,
`membraneThicknessForStraightening`, `bool_loadLabelmap` and the Cellpose environment are missing (issue I13).
The port saves the complete `settings.toml`, a `run.json` (version, counts, Cellpose version and device, timing) and
`log.txt`, plus `labelmap_final.tif` (cleaned labels) and `labelmap_classes.tif` (1 = immune cell, 2 = tumor
cell). In the saved `settings.toml` the channels given by name are replaced by their numbers.

Settings file ↔ macro parameters:

| Macro | Settings |
|---|---|
| `inputImageStack`, `output` | `input.image`, `input.output` (results in `<output>/<image stem>/`) |
| `bool_limitNrCells`, `limitNrCells` | `input.max_events` |
| `ch_BF`, `ch_MOI`, `ch_tumor`, `ch_Tcell` | `channels.brightfield`, `.marker_of_interest`, `.tumor`, `.immune` (a list; number or name) |
| (2.8: phenotypic, functional immune channels) | `channels.immune`, `channels.functional` |
| (2.8: `thresholdingMethod` automatic / manual) | `classification.threshold_method`, `classification.thresholds` |
| (2.8: `removeBadDetections`, constant) | `classification.artefact_min_positive` |
| (2.8: several input stacks, concatenated) | a project with `classification.pool = "experiment"` |
| `markerLoc`, `CellOfInterest` | `channels.marker_location`; `CellOfInterest` follows from the data (removed 2026-11-03) |
| `CellposeModel`, `cellDiameter` | `segmentation.model`, `segmentation.diameter` |
| `bool_loadLabelmap`, `labelmapFile` | `segmentation.method = "labelmap"`, `segmentation.labelmap_file` |
| `erosionSteps`, `minCellSize`, `maxCellSize` | `segmentation.erosion`, `.min_cell_size`, `.max_cell_size` |
| (new) | `segmentation.exclude_edge_cells` |
| `CellposeFlowThreshold` (constant) | `segmentation.flow_threshold` |
| `membraneThickness`, `membraneThicknessForStraightening` | `membrane.thickness`, `membrane.straightening_width` |
| `CellOfInterest` (the interface with the other class) | the cell type that carries the marker (from the data), or the pair `membrane.cells`, `membrane.partners` |
| (new) | `membrane.inside_only` |
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
