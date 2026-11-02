# Segmentation

The cells in every event are segmented with Cellpose. The resulting masks are then cleaned up, and the cells are
measured.

## The image for Cellpose

Cellpose does not segment each event on its own. Instead, it receives one **montage** of all events of a dataset,
placed side by side in a grid that is roughly square. What the montage contains depends on the setting
*Tumor / immune markers*:

- **nuclear** (the markers stain the nucleus or the cytoplasm): the montage has two channels. Channel 1 is the local
  variance of the brightfield image, which outlines the cells. Channel 2 is the sum of the normalized tumor and immune
  markers.
- **membrane**: the montage is a single image, the sum of the normalized markers and the normalized brightfield
  variance.

*Cell diameter (px)* is the typical cell diameter in pixels. *Model* is `cyto3` (Cellpose 3), `cpsam` (Cellpose 4), or
the path of a model of your own.

??? note "Technical: why a montage, and why a grid"

    A single ImageStream event shows a few cells on an empty background. That looks quite unlike the microscopy images
    on which the Cellpose models were trained. When the events were segmented one by one, the masks were worse and
    segmentation took longer. The Fiji macro put all events in a single row. A grid of about equal height and width is
    about 30 % faster and gives nearly the same masks (a median IoU of about 0.97 compared with the row, on validation
    data). Setting *Montage columns* to the number of events gives the single row of the macro.

    Normalization: as in the macro, each marker is divided by the median, over the events, of its Otsu threshold in
    each event. When there are several immune markers, each is normalized separately and the results are summed.

## Sizes in pixels

All sizes are given in pixels of the export. To convert a size in µm to pixels, divide it by the pixel size (1, 0.5
or 0.33 µm at 20×, 40× or 60×).

| Pixel size | 15 µm cell: diameter | its area | 25 µm cell: diameter | its area |
|---|---|---|---|---|
| 1 µm (20×) | 15 px | about 180 px | 25 px | about 490 px |
| 0.5 µm (40×) | 30 px | about 700 px | 50 px | about 2000 px |
| 0.33 µm (60×) | 45 px | about 1600 px | 75 px | about 4400 px |

- **Cell diameter (px)**: if the two cell types differ in size, start with a value between the two sizes, closer to
  the smaller one. Then check in the dashboard that large cells are not split and small cells are not merged.
  Alternatively, set it to 0 to let Cellpose 3 estimate the diameter.
- **Min. / Max. cell size (px)**: these are areas after shrinking. Keep the minimum well below the area of the smallest
  cells (about a third of it), and the maximum above the area of the largest cells.
- **Membrane / interface width (px)**: when left empty (the default), it is about 1.5 µm, which is 2 px at 20×, 3 px
  at 40× and 5 px at 60×. The same width is the distance within which two cells count as touching (see
  [Cell contacts](interactions.md)).

## Clean-up

1. A cell that Cellpose drew across the border between two events in the montage is split into one cell per event.
2. **Shrink cells (px)**: the cells are shrunk by this number of pixels where they border the background, but not
   where they border each other.
3. **Min. / Max. cell size (px)**: cells whose area after shrinking lies outside this range are removed.

??? note "Technical: details of the clean-up"

    Shrinking is a 4-connected binary erosion of all cells together, multiplied by the label image, as CLIJ2 does it
    in the macro. As a result, touching cells stay in contact, which the contacts and interfaces depend on. As in the
    macro, the size filter is only applied when the cells are shrunk (*Shrink cells* > 0). The labels are numbered
    1…N in the order of the events. Each cell's label is unique within the dataset and identifies the cell in the
    tables.

## Cells at the image edge

A cell that is cut by the edge of its event's own image (not of the padded frame) gets `Edge` = 1 in the cell table.
Its area and outline are incomplete, and so is its membrane. Tick *Leave out cells cut by the image edge* to remove
these cells before the classification.

## Measurements

Every cell is measured in every loaded channel. For each channel, the cell table holds the mean intensity minus the
channel's background (`Mean ch<N>`). The area of each cell is given in pixels and in µm². The normalized tumor and
immune markers give `Tumor intensity` and `Immune intensity`, which are used for the classification.

## Running again without Cellpose

Every run stores the Cellpose result of each dataset (`labelmap_cells.tif`) together with a fingerprint of what
produced it (`labelmap_cells.json`). When a dataset is run again and the fingerprint has not changed, Cellpose is
skipped, and the log says so. The clean-up, the measurements and all later steps are still done again. Changing a
threshold, the membrane pair or the clean-up settings (shrinking, cell sizes) therefore does not require a new
segmentation. Switch off *Reuse the previous Cellpose result if nothing changed* (`segmentation.reuse_labels`) to run
Cellpose every time.

With `segmentation.method = "labelmap"` in a settings file, the labels of another run (`labelmap_cells.tif`, or the
macro's `labelmap_cells.zip`) are loaded instead of running Cellpose.

??? note "Technical: what the fingerprint covers"

    The fingerprint is a hash of the image given to Cellpose and of the Cellpose settings. The image is the montage of
    brightfield and markers, so the hash covers the data, the channels, the number of events, the marker location and
    the number of montage columns. The settings are the model, the diameter, the flow and cell probability thresholds,
    and the Python that runs Cellpose. The GPU setting is not part of it. Results of runs made before the
    fingerprint was introduced have no fingerprint, so Cellpose runs once more the next time.
