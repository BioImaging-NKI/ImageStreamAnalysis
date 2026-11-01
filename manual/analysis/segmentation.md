# Segmentation

The cells of every event are segmented with Cellpose, then cleaned up and measured.

## The image for Cellpose

Cellpose gets one **montage** of all events of a dataset, side by side in an about-square grid, not every event on its
own. Its input depends on *Tumor / immune markers*:

- **nuclear** (the markers stain the nucleus or the cytoplasm): channel 1 is the local variance of the brightfield
  (which outlines the cells), channel 2 the sum of the normalized tumor and immune markers.
- **membrane**: one image, the sum of the normalized markers and of the normalized brightfield variance.

*Cell diameter (px)* is the typical diameter in pixels; *Model* is `cyto3` (Cellpose 3) or `cpsam` (Cellpose 4), or
the path of a model of your own.

??? note "Technical: why a montage, and why a grid"

    A single ImageStream event, a few cells on an empty background, looks unlike the microscopy images that the
    Cellpose models are trained on; segmented one by one, the events gave worse masks and took longer. The Fiji macro
    put all events in one row; a grid of about the same height and width is about 30 % faster and gives nearly the
    same masks (median IoU about 0.97 with the row on validation data). *Montage columns* = the number of events
    gives the macro's row.

    Normalization: each marker is divided by the median, over the events, of its Otsu threshold per event (as the
    macro). Several immune markers are each normalized and then summed.

## Sizes in pixels

All sizes are in pixels of the export: a size in µm divided by the pixel size (1, 0.5 or 0.33 µm at 20×, 40× or 60×).

| Pixel size | 15 µm cell: diameter | its area | 25 µm cell: diameter | its area |
|---|---|---|---|---|
| 1 µm (20×) | 15 px | about 180 px | 25 px | about 490 px |
| 0.5 µm (40×) | 30 px | about 700 px | 50 px | about 2000 px |
| 0.33 µm (60×) | 45 px | about 1600 px | 75 px | about 4400 px |

- **Cell diameter (px)**: with two cell types of different size, start between the two, closer to the smaller one,
  and check in the dashboard that large cells are not split and small ones not merged (or let Cellpose 3 estimate
  it: 0).
- **Min. / Max. cell size (px)**: areas after shrinking; keep the minimum well below the area of the smallest cells
  (about a third of it) and the maximum above the largest.
- **Membrane / interface width (px)**: empty (the default) is about 1.5 µm: 2 px at 20×, 3 px at 40×, 5 px at 60×.
  It is also the distance at which two cells count as touching (see [Cell contacts](interactions.md)).

## Clean-up

1. A cell that Cellpose drew across the border of two events in the montage is split into one cell per event.
2. **Shrink cells (px)**: the cells are shrunk from the background by this many pixels (not from each other).
3. **Min. / Max. cell size (px)**: cells outside the range, after shrinking, are removed.

??? note "Technical: details of the clean-up"

    The shrinking is a 4-connected binary erosion of all cells together, multiplied with the labels, as CLIJ2 does in
    the macro: touching cells stay touching, which the contacts and interfaces rely on. As in the macro, the size
    filter is only applied when the cells are shrunk (*Shrink cells* > 0). Labels are numbered 1…N in the order of
    the events; a cell's label is unique in the dataset and names it in the tables.

## Cells at the image edge

A cell cut by the edge of its event's own image (not of the padded frame) gets `Edge` = 1 in the cell table. *Leave out
cells cut by the image edge* removes them before the classification. Their area and outline are incomplete, and so
is their membrane.

## Measurements

Every cell is measured in every loaded channel: the mean minus the channel's background (`Mean ch<N>`), and its area
in pixels and µm². The normalized tumor and immune markers give `Tumor intensity` and `Immune intensity`, used for the
classification.

## Running again without Cellpose

Every run keeps the Cellpose result of a dataset (`labelmap_cells.tif`) with a fingerprint of what produced it
(`labelmap_cells.json`). The next run of the dataset skips Cellpose when the fingerprint is the same, and says so in
the log; clean-up, measurements and everything after them are done again. So changing a threshold, the membrane pair
or the clean-up settings (shrinking, cell sizes) costs no new segmentation. *Reuse the previous Cellpose result if
nothing changed* (`segmentation.reuse_labels`) switched off always runs Cellpose.

With `segmentation.method = "labelmap"` in a settings file, the labels of another run (`labelmap_cells.tif`, or the
macro's `labelmap_cells.zip`) are loaded instead of running Cellpose.

??? note "Technical: what the fingerprint covers"

    The fingerprint is a hash of the image given to Cellpose (the montage of brightfield and markers, so it covers the
    data, the channels, the number of events, the marker location and the montage columns) and of the model, the
    diameter, the flow and cell probability thresholds and the Python of Cellpose. The GPU setting is not part of it.
    Results of runs before this was added have no fingerprint, so their next run uses Cellpose once more.
