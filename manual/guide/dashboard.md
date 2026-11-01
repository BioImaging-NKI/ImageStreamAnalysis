# The dashboard

Every analysed dataset gets a `dashboard.html` in its folder. It opens from the overview of the import page, or
directly from disk (double-click): it needs no server and no internet, so a folder can be sent to a colleague. The
event images are in `dashboard_files/` next to it.

Everything is **linked**: hovering a cell in the events, a plot or a table highlights it everywhere; clicking selects
it and shows it in the *Cell* card. *Esc* clears the selection. `dashboard.html#cell=12` or `#event=5` in the address
selects that cell (row of the cell table) or event when the page opens.

![The dashboard](../images/dash_overview.png)

## Top: counts and filters

The header counts events, cells, cell types, events with one tumor and one immune cell, touching pairs, artefacts and
cells with an interface. The filter row below applies to every view:

- **Events**: all, those with immune and tumor cells, one tumor + one immune cell, those with an interface, those with
  touching cells. A contact chosen in the [interaction matrix](#interactions) adds a chip here (✕ clears it).
- **Sort**: the order of the events.
- **Size**: the size of the event images.

The **?** button behind a heading, and the **? Help** button at the top, open the part of this manual about it in a
panel at the side; <kbd>F1</kbd> or <kbd>?</kbd> does the same for the card under the mouse. The panel works offline.

![The help of the subtype markers, next to the card](../images/dash_help.png){ width="600" }

## Events

The events as a gallery, with the cells outlined in the colour of their type and the interfaces in aqua. Only the
visible events are drawn, so large datasets stay fast.

## Channels and overlays

The column left of the gallery switches what is drawn in the event images:

- **Channels**: one row per channel, with a tick box and its brightness. Double-click a channel's name to show only
  that channel. Brightfield starts at half brightness, under the fluorescence.
- **Overlays**: the **outlines** of the cells (in the colour of their type), the **interfaces** and the **contacts**:
  a white line, with a dark edge, between the centres of every two cells that touch.

## Cell

The selected cell: its event, measurements, subtype markers (value, positive or negative), the cells it touches, and,
with the membrane analysis, its membrane profile aligned on the interface.

![A selected cell](../images/dash_cell.png){ width="420" }

## Subtype markers

Shown when the panel has subtype markers. Per marker, a histogram of the immune cells (on log(intensity − min + 1)),
the threshold of the run (dashed) and a slider. **Dragging the slider changes the subtypes everywhere on the page** at
once: the colours, the plots, the contacts, the membrane pairs and the tables. The count next to the marker says how
many immune cells are positive.

![Subtype markers](../images/dash_markers.png){ width="420" }

*Copy the thresholds as settings* copies a `[classification]` line with the thresholds to the clipboard (and shows
it): paste it into the project or settings file and run again to make them the thresholds of the tables. *Reset all*
returns to the run's thresholds.

??? note "Technical: what the slider changes, and what not"

    The slider re-applies the threshold to the immune cells in the page. It does not repeat the tumor / immune
    clustering, and the artefacts stay as computed: both would need the whole classification again (the artefact rule
    uses thresholds over *all* detections, see [Cell types](../analysis/classification.md#artefacts)). In the run, a
    threshold set by hand (`classification.thresholds`) is used for both: for the subtypes and for the artefact rule.

    To run again without segmenting again, set `segmentation.method = "labelmap"` and `segmentation.labelmap_file` to
    the dataset's `labelmap_cells.tif` (one dataset), or keep the project as it is (Cellpose runs again; the result is
    the same for the same settings).

## Cell classes

Area against the marker ratio (immune / tumor marker), log–log: the space in which the k-means separates tumor cells
from immune cells. Clicking a type in the legend brings it forward.

![Cell classes](../images/dash_classes.png){ width="420" }

## Marker of interest at the interface

Per event, the mean of the marker of interest at the interface against the rest of the membrane, for the cell type
that expresses it (the per-event table of the Fiji macro). Points above the dashed line: more at the interface.

## Interactions

![Interactions](../images/dash_interactions.png)

- **Partner preference** (left): rows are cell types, columns the type they touch. The colour and the number:
  contacts observed ÷ expected by chance, within the events of three or more cells (blue: fewer than chance, red:
  more). *fraction touching* shows instead the share of the row's cells that touch one of the column's type.
  Hovering a square gives the counts. **Clicking a square shows only the events with that contact**, with the
  contacts drawn.
- **Share of contacts** (right): per cell type, whom it touches, as a share of its contacts or as numbers; on the right
  the share of its cells that touch any cell.

How contacts and chance are defined, and what they can and cannot tell, is in
[Cell contacts](../analysis/interactions.md).

![The events with a contact of tumor cells and Subtype A cells](../images/dash_contact_filter.png){ width="560" }

## Membrane enrichment

Shown when the membrane analysis was run. **Membrane of** and **interface with** choose the pair: a cell type, an
immune subtype, or all immune cells positive for a marker. The page computes the enrichment for that pair at once.

![Membrane enrichment](../images/dash_membrane.png)

- **Histogram**: the enrichment at the interface (aqua) of every cell with an interface, against random control
  windows on cells without an interface (dashed), in standard deviations above the mean of the free membrane.
- **Average profile**: the profiles of the cells with an interface, centred on the interface.
- **Heatmap**: one row per cell with an interface, the strongest on top.
- **The badge** next to the menus: whether the marker of interest is on the measured cells or also on their partners.
  A red badge means the interface may show the partner's membrane. Read
  [Membrane enrichment](../analysis/membrane.md#where-the-marker-is-expressed) before trusting a result with an orange
  or red badge.

??? note "Technical: the page and the tables"

    The page computes the pair from every cell's membrane profile and neighbour track, as the program does for the
    tables. Only the random control windows differ (the page uses a fixed seed, the run an unseeded generator); their
    summaries agree closely. The tables (`Membrane intensity analysis.tsv`) hold the pair of the settings.

## Tables

*Cells*, *Events* and *Contacts*, sortable by clicking a column header. Rows follow the filters; at most 2000 are
shown (sort to see others). The full tables are the TSV files in the [output folder](output-folder.md).
