# The output folder

## The project folder

| File | Contents |
|---|---|
| `project.toml` | the project: the exports, channel names, roles and settings (see [Settings files](../reference/command-line.md)) |
| `overview.tsv` | one row per dataset (the [overview](overview.md)) |
| `interactions.tsv` | the contact tables of all datasets, and their sums per experiment |
| `membrane pairs.tsv` | the membrane pairs of all datasets |
| `condition comparison.tsv` | the comparison of the conditions (only when every dataset has a condition and a replicate): one row per measure, pair and condition, with the change, its 95 % interval, *p*, *q*, the value per replicate, and notes; see [Comparing conditions](../analysis/comparing-conditions.md) |
| `subset preference.tsv` | the subset preference of all datasets (with intervals), and its sums per experiment (without intervals) |
| `<experiment>/<dataset>/` | the results of one dataset (see below) |

## The folder of a dataset

The tables are tab-separated text files (`.tsv`) with 3 decimals. `NaN` marks a value that does not exist, such as the
mean of no values. For an IDEAS export, the tables with one row per cell or per event end with the column `Object ID`,
which holds the IDEAS object number.

| File | Contents |
|---|---|
| `Cell classification.tsv` | one row per cell: the measurements, cell type, subtype and neighbors (see below) |
| `All cell combinations.tsv` | one row per event: the number of cells of each type, the marker of interest at the interfaces (left out when there is no marker of interest), whether a tumor cell and an immune cell touch (`Tumor-immune interface`), and the contacts |
| `Single Tumor Single Immune cell.tsv` | the events with exactly one tumor cell and one immune cell that touch |
| `Marker thresholds.tsv` | the threshold of every subtype marker, both for the artefact rule and on the immune cells, as background-subtracted intensity (on the same scale as `Mean ch<N>`), and how it was found (`how`) |
| `Interactions.tsv` | the contacts per cell type and neighbor type, compared with chance |
| `Membrane pairs.tsv` | every pair of *Membrane of* and *Interface with*: the cells with an interface, the median enrichment and ratio, and the marker check (only with the membrane analysis) |
| `Subset preference.tsv` | per immune subtype and partner: the contacts compared with those when the subtypes are shuffled among the immune cells of each event, with the 95 % interval (only with two or more immune subtypes) |
| `Multiplets.tsv` | the number of events with n cells and with n contacts |
| `Membrane intensity analysis.tsv` | one row per measured cell: the interface and membrane means, the enrichment, and the control |
| `Average normalized membrane profile values - interfaces centered at 20.tsv` | the average aligned profile |
| `Histogram_interface_stdDevs_above_mean.tsv`, `Histogram_control_stdDevs_above_mean.tsv` | the histograms |
| `Control windows stdDevs above mean.tsv` | every random control window |
| `membrane_profiles_normalized.tif` | the aligned profiles, with one row for each row of the membrane table |
| `labelmap_cells.tif`, `labelmap_cells.json` | the labels as Cellpose made them (the montage), and what produced them; the next run reuses them when nothing has changed (see [Segmentation](../analysis/segmentation.md#running-again-without-cellpose)) |
| `labelmap_final.tif`, `labelmap_classes.tif` | the cleaned-up labels per event; the cell type per pixel (1 immune, 2 tumor, 3 artefact) |
| `dashboard.html`, `dashboard_files/` | the [dashboard](dashboard.md) |
| `settings.toml` | the complete settings of the run (with channels by number) |
| `run.json`, `summary.json`, `log.txt` | the version, counts, Cellpose version and device, timing and marker check; the log |

The membrane files exist only when the membrane analysis was run.

## The cell table

| Column | Meaning |
|---|---|
| `Immune intensity`, `Tumor intensity` | the mean of the normalized immune markers (summed) and of the normalized tumor marker in the cell |
| `Marker ratio`, `Marker ratio over Area` | the immune intensity divided by the tumor intensity; the same ratio divided by the area (the k-means uses these) |
| `Area`, `Area (µm²)` | the area in pixels; the area in square micrometers |
| `Class`, `Cell type` | 0 tumor, 1 immune, 2 artefact; Tumor, Immune, Artefact |
| `Frame`, `Label` | the event (1-based) and the cell's label (unique within the dataset) |
| `Edge` | 1 when the cell is cut by the edge of its event image |
| `Mean ch<N> (<name>)` | the mean of every loaded channel, minus the background of that channel |
| `<marker>+` | 0 or 1 for each subtype marker (using the threshold of the immune cells) |
| `Positive markers` | the number of subtype markers in which the detection is positive (using the thresholds of the artefact rule) |
| `Subtype` | the immune subtype (for example "A+ B+" or "negative"); empty for other cells |
| `Neighbours`, `Touching (labels)`, `Neighbour types` | the cells it touches: how many, which ones (by label), and their types |
| `Object ID` | the IDEAS object number of the event |

??? note "Technical: names that differ from the Fiji macro"

    Where the macro v1.8 says *T cell*, this program says *immune cell* (as ClassifyCellTypes 2.8 does), because the
    immune side can be defined by several markers. The macro's table `Single Tumor Single T cell` is therefore called
    `Single Tumor Single Immune cell` here, and its columns `T-cell intensity`, `T cells`, `mean Tcell interface` and
    `T cell membrane ratio` are called `Immune intensity`, `Immune cells`, `mean immune interface` and
    `Immune membrane ratio`. The other columns keep the macro's names, and new columns come after them.
