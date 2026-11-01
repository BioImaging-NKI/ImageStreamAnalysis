# The output folder

## The project folder

| File | Contents |
|---|---|
| `project.toml` | the project: exports, channel names, roles and settings (see [Settings files](../reference/command-line.md)) |
| `overview.tsv` | one row per dataset (the [overview](overview.md)) |
| `interactions.tsv` | the contact tables of all datasets, and summed per experiment |
| `subset preference.tsv` | the subset preference of all datasets (with intervals), and summed per experiment (without) |
| `<experiment>/<dataset>/` | the results of one dataset (below) |

## The folder of a dataset

Tables are tab-separated text (`.tsv`) with 3 decimals; `NaN` marks a value that does not exist (a mean of nothing).
For an IDEAS export, the tables per cell and per event end with the column `Object ID`, the IDEAS object number.

| File | Contents |
|---|---|
| `Cell classification.tsv` | one row per cell: measurements, cell type, subtype, neighbours (below) |
| `All cell combinations.tsv` | one row per event: cells per type, marker of interest at the interfaces (left out without a marker of interest), whether a tumor and an immune cell touch (`Tumor-immune interface`), contacts |
| `Single Tumor Single Immune cell.tsv` | the events with exactly one tumor and one immune cell that touch |
| `Marker thresholds.tsv` | the threshold of every subtype marker, on all detections and on the immune cells, in background-subtracted intensity (as `Mean ch<N>`), and how it was found (`how`) |
| `Interactions.tsv` | contacts per cell type and neighbour type, with the comparison with chance |
| `Subset preference.tsv` | per immune subtype and partner: contacts against the subtypes shuffled among the immune cells of each event, with the 95 % interval (with two or more immune subtypes) |
| `Multiplets.tsv` | the number of events with n cells and with n contacts |
| `Membrane intensity analysis.tsv` | one row per measured cell: interface and membrane means, enrichment, control |
| `Average normalized membrane profile values - interfaces centered at 20.tsv` | the average aligned profile |
| `Histogram_interface_stdDevs_above_mean.tsv`, `Histogram_control_stdDevs_above_mean.tsv` | the histograms |
| `Control windows stdDevs above mean.tsv` | every random control window |
| `membrane_profiles_normalized.tif` | the aligned profiles, one row per row of the membrane table |
| `labelmap_cells.tif`, `labelmap_cells.json` | the labels as Cellpose made them (the montage), and what produced them: reused by the next run when nothing changed (see [Segmentation](../analysis/segmentation.md#running-again-without-cellpose)) |
| `labelmap_final.tif`, `labelmap_classes.tif` | the cleaned labels per event; the cell type per pixel (1 immune, 2 tumor, 3 artefact) |
| `dashboard.html`, `dashboard_files/` | the [dashboard](dashboard.md) |
| `settings.toml` | the complete settings of the run (channels by number) |
| `run.json`, `summary.json`, `log.txt` | version, counts, Cellpose version and device, timing, the marker check; the log |

The membrane files exist only when the membrane analysis was run.

## The cell table

| Column | Meaning |
|---|---|
| `Immune intensity`, `Tumor intensity` | mean of the normalized immune markers (summed) and tumor marker in the cell |
| `Marker ratio`, `Marker ratio over Area` | immune / tumor intensity; the same divided by the area (the k-means uses these) |
| `Area`, `Area (µm²)` | pixels; square micrometres |
| `Class`, `Cell type` | 0 tumor, 1 immune, 2 artefact; Tumor, Immune, Artefact |
| `Frame`, `Label` | the event (1-based) and the cell's label (unique in the dataset) |
| `Edge` | 1: the cell is cut by the edge of its event image |
| `Mean ch<N> (<name>)` | the mean of every loaded channel, minus the channel's background |
| `<marker>+` | 0/1 per subtype marker (with the threshold of the immune cells) |
| `Positive markers` | the number of subtype markers the detection is positive in (thresholds over all detections: the artefact rule) |
| `Subtype` | the immune subtype ("A+ B+", "negative"); empty for other cells |
| `Neighbours`, `Touching (labels)`, `Neighbour types` | the cells it touches: how many, which (labels), and their types |
| `Object ID` | the IDEAS object number of the event |

??? note "Technical: names that differ from the Fiji macro"

    The macro v1.8 says *T cell* where this program says *immune cell* (as ClassifyCellTypes 2.8), because the
    immune side can be defined by several markers. Its table `Single Tumor Single T cell` is therefore
    `Single Tumor Single Immune cell` here, and its columns `T-cell intensity`, `T cells`, `mean Tcell interface`
    and `T cell membrane ratio` are `Immune intensity`, `Immune cells`, `mean immune interface` and
    `Immune membrane ratio`. The other columns keep the macro's names; new columns come after them.
