# The import page

The import page sets up a **project**. A project consists of one or more IDEAS exports, the names and roles of their
channels, and the analysis settings. From this page you also run the analysis and view the overview of the results.
Start it with `uv run imagestream` or `uv run imagestream-cli gui`. The page opens in your browser at once, with a
spinner while the program loads its analysis modules (up to about 20 s on a first start). It is served by
the program on your own computer, so no data leave your computer.

The page remembers the last saved project. *Save settings*, *Run* and *Quit* all save the project. If you reload the page
before saving, the browser restores the unsaved setup and shows a message about it. If you start the program a second
time, it opens the page that is already running. **Quit** (top right) saves the project and stops the program. If you
don't quit, the program stops by itself 15 minutes after the page was closed, but never during a run.

The **?** button next to a step or a group of settings, and the **Help** button at the top, open the matching part
of this manual in a panel at the side. Pressing <kbd>F1</kbd> or <kbd>?</kbd> does the same for the step under the
mouse pointer.

## 1. Project

Choose the output folder of the project here, or click *Browse…* to choose a folder. The results of each dataset go
to `<project>/<experiment>/<dataset>/`.

The top bar has *Load settings…*, which opens an existing `project.toml`, and *Save settings*, which writes the
project file (`project.toml` in the output folder) without running the analysis. A second *Save settings* button
sits at the bottom of the page, so that you don't have to scroll back up after changing the settings.

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
    | **MOI** | a marker of interest, which is measured along the membrane | one or more, or none |

    A channel can have two roles, for example a subtype marker that is also a marker of interest. With several markers
    of interest, the membrane enrichment is measured for each of them in the same way (see
    [Several markers of interest](../analysis/membrane.md#several-markers-of-interest)). Without a marker of interest,
    the analysis gives only the cell types, subtypes and contacts. In that case, leave *MOI* unticked and switch
    *Membrane profiles and enrichment* off.

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
- **Python with Cellpose**, **Model** and **Cell diameter (px)** control the segmentation. *Python with Cellpose* can
  be the `python.exe` of an environment or the folder of that environment (a venv, conda or pixi environment, or a
  project folder with `.venv` or `.pixi` in it). A chip in the heading of *Segmentation (Cellpose)* shows what the
  program found there: in green the Cellpose version and the name of the GPU, or a warning when Cellpose is missing
  or no GPU was found (Cellpose then runs on the processor, which is much slower). The check runs in the
  background each time the page opens. Sizes are given in pixels;
  see [Segmentation](../analysis/segmentation.md#sizes-in-pixels) for how to convert them from µm. **Reuse the previous
  Cellpose result if nothing changed** (on by default) skips Cellpose for a dataset when its Cellpose input and
  settings are the same as in its previous run.
- **Classify together**, **Threshold of the subtype markers** and **Artefact** control how the
  [cell types](../analysis/classification.md) are assigned. **Thresholds set by hand** lists the thresholds that were
  applied from a [dashboard](dashboard.md#subtype-markers); *Automatic again* removes them.
- **Membrane profiles and enrichment** and the membrane pairs control the
  [membrane enrichment](../analysis/membrane.md). Every marker of interest has a row of its own, for example
  *Marker: membrane of* … *interface with* …, so that each marker can be measured on its own pair. The rows follow the
  *MOI* ticks in step 2: a marker that you tick gets a row, and a marker that you untick loses its row and its pair.
  A renamed channel keeps its pair. Without a marker of interest there are no rows. *Membrane of* and *Interface with*
  are menus. Their first choice is automatic: *where the marker is expressed* means the cell type that carries the marker of interest, which the
  program finds in the data for each marker, and *the other cell type* means the other one. The other choices are the
  tumor cells, all immune cells, every combination of the subtype markers ticked in step 2 (for example *immune A+
  cells* or *immune negative cells*), and with two or more subtype markers each marker on its own (*immune cells A+
  (any other markers)*). These are the same choices as in the dashboard's membrane card, which shows only the
  combinations that occur in the data. The menus follow the roles: when you tick, untick or rename a subtype marker,
  they change with it. A pair that you choose applies to that marker only, and each side on its own: if you choose
  only *Membrane of*, *Interface with* stays automatic. The run check asks for at least one *MOI* channel when
  *Membrane profiles and enrichment* is on.

    A value that the menus do not offer, for example from an older project or a settings file edited by hand, is
    kept and shown as its own choice, marked *(from the settings file)*, so that saving the project does not change
    it. If it is not a selection of the panel, the list above the *Run* button says so, with the likely meaning: for
    example, `Immune A` is no subset, and the program asks whether you meant `A+` (every A-positive immune cell) or
    `Immune A+` (A positive only). The run does not start until you choose a pair in the menu. A valid subset that
    has no cells in a dataset only gives a warning in the log, and that dataset gets no membrane files for it.

    In `project.toml` the pairs are tables by marker name, for example `[membrane.pairs.'Marker']` with
    `cells = ['Tumor']` and `partners = ['Immune Subtype A+']`; the menus write one entry per side, or none for the
    automatic choice. An older project with one pair for all markers (`membrane.cells` and `membrane.partners`) is
    shown with that pair in the row of every marker, and saved as one entry per marker.

[All settings](../reference/settings.md) lists every setting with its default value.

### Pools

The setting *Classify together (one set of thresholds)* has two options. With **experiment** (the default), the cells
of all datasets with the same experiment label are classified together. As a result, "positive for a marker" means
the same in every dataset of that experiment. With **dataset**, each dataset is classified on its own. See
[Cell types](../analysis/classification.md#several-datasets-pools).

## 4. Run

The list above *Run all datasets* names anything that is still missing, such as an output folder, a pixel size, or a
channel for a role. For the immune-cell markers, one ticked *Immune* channel in each export is enough: the markers are
summed, so an export may lack a marker that another export has. The button works only when this list is empty. The
run first saves the project. The log shows the progress of the run, and *Stop* ends it. If a dataset fails, the error
is logged and the other datasets are still analyzed.

To run **one dataset** (again), click *Run* in its row of the dataset table in step 2, or *run again* in its row of the
overview. The other datasets keep their results, and the overview is written again from all of them. With *Classify
together* set to *experiment*, a dataset is classified together with the other datasets of its experiment. Its whole
experiment then runs, so that all of them keep one set of thresholds; the tooltip of the button says so. Cellpose is
skipped for every dataset whose Cellpose input and settings have not changed.

## 5. Overview

The overview shows one row per analyzed dataset, with links to its [dashboard](dashboard.md) and its folder. Below
the table, it shows the partner preference of each experiment and, if the datasets have conditions, the comparison of
the conditions. See [The project overview](overview.md).
