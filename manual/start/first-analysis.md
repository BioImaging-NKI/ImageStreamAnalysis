# Your first analysis

This walk-through analyses the example data: a synthetic IDEAS export with two datasets ("Donor 1" and "Donor 2") of
tumor cells and immune cells, with two immune subtype markers and a marker of interest on the tumor cells. The data are
made by a script, so that the true cell types are known; they are not measurements.

## 1. Get the example

Download `example_export.zip` from the
[example-data release](https://github.com/BioImaging-NKI/ImageStreamAnalysis/releases/tag/example-data) and put it in
a folder of your own, for example `D:\data`. It is an export as IDEAS writes it: one folder per dataset, one TIFF per
object and channel (`<object>_Ch<N>.ome.tif`).

## 2. Open the import page

```powershell
uv run imagestream
```

The import page opens in your browser. It has five steps, from top to bottom.

**Step 1, Project.** Choose an empty folder for the results, for example `D:\data\example project`. The project file
(`project.toml`) and all results go there.

## 3. Add the export and name the channels

**Step 2, IDEAS exports.** Click *Add zip files…* and choose `example_export.zip`. The page reads the zip and shows
its datasets and channels, with a few example objects of each channel to recognise them by.

![The IDEAS exports step with the example](../images/gui_exports.png)

1. **Pixel size**: *40× (0.5 µm)*. IDEAS files do not contain it.
2. **Datasets**: both are ticked. The *Experiment* label groups datasets: datasets with the same label are classified
   together and shown together in the overview. Both are *Example* here.
3. **Channels**: type a name for every channel that is needed, and tick its role:

    | Channel | Name | Role |
    |---|---|---|
    | Ch1 | BF | BF (brightfield) |
    | Ch2 | Marker | MOI (marker of interest) |
    | Ch3 | Tumor | Tumor |
    | Ch4 | Immune | Immune |
    | Ch7 | Subtype A | Subtype |
    | Ch8 | Subtype B | Subtype |

    Ch9 is not needed: leave it without a name, and it is not loaded.

Roles go with the *name*: a second export of the same panel, even with its channels numbered differently, only needs
the same names. [The import page](../guide/import-page.md) describes every step.

## 4. Settings

**Step 3, Analysis settings.** For the example, change:

- *Marker of interest expressed in*: **Tumor cell** (the marker is on the tumor cells' membrane);
- *Cell diameter (px)*: **20**;
- *Membrane profiles and enrichment*: **on**.

Leave the rest. If Cellpose is installed in another environment, fill in *Python with Cellpose* (see
[Installation](install.md#3-cellpose)).

![The analysis settings](../images/gui_settings.png)

## 5. Run

**Step 4, Run.** The list above the button says what is still missing, if anything. Click *Run all datasets*. The log
shows the progress: first every dataset is segmented and measured, then the datasets of each experiment are
classified together, then every dataset is finished. The example takes about a minute with a GPU.

## 6. The overview

**Step 5, Overview.** One row per dataset: counts of cells and cell types, contacts, and the membrane results, with a
link to each dataset's dashboard. Below it, the partner preference of each experiment.

![The project overview](../images/gui_overview.png)

## 7. The dashboard

Click *dashboard* in the row of Donor 1. The [dashboard](../guide/dashboard.md) shows every event with the cells
outlined in the colour of their type, and the results linked: hover or click a cell anywhere and it is highlighted in
every view.

![The dashboard of Donor 1](../images/dash_overview.png)

Things to try:

- **Subtype markers**: drag the threshold of *Subtype A+* and watch the subtypes change everywhere.
- **Interactions**: click a square of the matrix, for example *tumor* × *Subtype A+*, to see only the events with
  that contact; tick *contacts* to draw them.
- **Membrane enrichment**: choose *interface with* **immune cells Subtype A+**: in the example, the marker of interest
  gathers at the contacts with those cells, not with the others.

The tables and files that the run wrote are described in [The output folder](../guide/output-folder.md).
