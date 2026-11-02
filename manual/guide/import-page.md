# The import page

The import page sets up a **project**. A project consists of one or more IDEAS exports, the names and roles of their
channels, and the analysis settings. From this page you also run the analysis and view the overview of the results.
Start it with `uv run imagestream` or `uv run imagestream-cli gui`. The page opens in your browser, but it is served by
the program on your own computer, so no data leave your computer.

The page remembers the last saved project. *Save*, *Run* and *Quit* all save the project. If you reload the page
before saving, the browser restores the unsaved setup and shows a message about it. If you start the program a second
time, it opens the page that is already running. **Quit** (top right) saves the project and stops the program. If you
don't quit, the program stops by itself 15 minutes after the page was closed, but never during a run.

The **?** button next to a step or a group of settings, and the **? Help** button at the top, open the matching part
of this manual in a panel at the side. Pressing <kbd>F1</kbd> or <kbd>?</kbd> does the same for the step under the
mouse pointer.

## 1. Project

Choose the output folder of the project here. Click *Browse…* to choose a folder, or *Open project…* to open an
existing `project.toml`. *Save* writes the project file without running the analysis. The results of each dataset go
to `<project>/<experiment>/<dataset>/`.

## 2. IDEAS exports

Click *Add zip files…* to add exports, or paste the path of a zip file and click *Add*. For each export, set the
following:

- **Pixel size**: choose 20×, 40× or 60× (1, 0.5 or 0.33 µm per pixel), or *other…* to type a value. The IDEAS files
  don't contain the pixel size, so you have to set it. These are the nominal values; check them for your instrument.
- **Datasets**: each first-level folder in the zip is one dataset, listed in a table. Tick *Use* to include a dataset
  in the analysis. *Objects* shows how many objects the dataset contains, and *Complete* how many of them have an
  image for every channel that is needed; incomplete objects are left out. *Experiment* groups datasets that are
  classified together (see [Pools](#pools) below); it starts as the name of the zip. *Condition* (for example control
  or treated) and *Replicate* (the donor or independent experiment the dataset comes from) are needed only if you want
  to compare conditions; see [Comparing conditions](../analysis/comparing-conditions.md#labeling-the-datasets).
- **Channels**: the example images show a few objects in each channel. Give every channel that is needed a **name**,
  and tick its **role**:

    | Role | Meaning | How many |
    |---|---|---|
    | **BF** | brightfield | one |
    | **Tumor** | the tumor-cell marker | one |
    | **Immune** | an immune-cell marker; if there are several, they are summed | one or more |
    | **Subtype** | an immune subtype marker; each cell is called positive or negative for it | none or more |
    | **MOI** | the marker of interest, which is measured along the membrane | one, or none |

    A channel can have two roles, for example a subtype marker that is also the marker of interest. Without a marker
    of interest, the analysis gives only the cell types, subtypes and contacts. In that case, leave *MOI* unticked and
    switch *Membrane profiles and enrichment* off.

    When two or more exports have been added, a *Copy setup to all* button appears at the bottom of each export's
    card, below the channel table. It copies the pixel size and the channel names of that export to the other
    exports. The roles follow the names (see below), so they need no copying.

![Channels with their names and roles](../images/gui_channels.png)

*Remove*, in the title bar of an export, removes the export from the project, but not from the disk.

??? note "Technical: roles go with the name"

    The roles are stored once for the whole project, by marker name (`[channels]` in `project.toml`). Each export only
    maps its channel numbers to names. As a result, exports of the same panel whose channels are numbered differently
    (because of another instrument setup or another IDEAS template) need no separate roles, only the right names. A
    role tick box is grayed out until the channel has a name. When you rename a channel, its roles move to the new
    name, unless another export still uses the old name.

    Projects written before November 2026 stored the roles per export, by channel number. When such a project is
    opened, the roles are converted to roles by name if all exports agree. Otherwise, each export keeps its own roles
    by number, which are shown under its channel table.

## 3. Analysis settings

These settings apply to all datasets of the project. The most important ones are:

- **Tumor / immune markers**: *nuclear* or *membrane*. This tells the program where the tumor and immune markers
  stain, which determines how the image for Cellpose is made (see [Segmentation](../analysis/segmentation.md)).
- **Marker of interest expressed in**: the cell type that carries the marker of interest. The membrane enrichment is
  measured on these cells (see [Membrane enrichment](../analysis/membrane.md)).
- **Python with Cellpose**, **Model** and **Cell diameter (px)** control the segmentation. Sizes are given in pixels;
  see [Segmentation](../analysis/segmentation.md#sizes-in-pixels) for how to convert them from µm. **Reuse the previous
  Cellpose result if nothing changed** (on by default) skips Cellpose for a dataset when its Cellpose input and
  settings are the same as in its previous run.
- **Classify together**, **Threshold of the subtype markers** and **Artefact** control how the
  [cell types](../analysis/classification.md) are assigned. **Thresholds set by hand** lists the thresholds that were
  applied from a [dashboard](dashboard.md#subtype-markers); *Automatic again* removes them.
- **Membrane profiles and enrichment**, **Membrane of** and **Interface with** control the
  [membrane enrichment](../analysis/membrane.md).

[All settings](../reference/settings.md) lists every setting with its default value.

### Pools

The setting *Classify together (one set of thresholds)* has two options. With **experiment** (the default), the cells
of all datasets with the same experiment label are classified together. As a result, "positive for a marker" means
the same in every dataset of that experiment. With **dataset**, each dataset is classified on its own. See
[Cell types](../analysis/classification.md#several-datasets-pools).

## 4. Run

The list above *Run all datasets* names anything that is still missing, such as an output folder, a pixel size, or a
channel for a role. The button works only when this list is empty. The run first saves the project. The log shows the
progress of the run, and *Stop* ends it. If a dataset fails, the error is logged and the other datasets are still
analyzed.

## 5. Overview

The overview shows one row per analyzed dataset, with links to its [dashboard](dashboard.md) and its folder. Below
the table, it shows the partner preference of each experiment and, if the datasets have conditions, the comparison of
the conditions. See [The project overview](overview.md).
