# From the Fiji macros

ImageStream Analysis for Python is a port of two Fiji macros, which are kept in the
[`legacy/`](https://github.com/BioImaging-NKI/ImageStreamAnalysis/tree/main/legacy) folder of the repository:

- **`Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm`** separates tumor cells and T cells
  by k-means, measures the interface in each event, and analyzes the membrane enrichment.
- **`ClassifyCellTypes_2_8.ijm`** uses several immune markers, assigns immune subtypes by thresholds, finds the cells
  that touch, and has a class for detections that are not cells.

The port reads the IDEAS exports directly, so it replaces ImageStreamCombiner. It also combines the two macros: the
analysis of macro 1.8 is the special case of one immune marker and no subtype markers. The port was validated against
the output of macro 1.8. The cells and the per-event tables match to three decimals, with 100 % agreement on the
classes, and the results of the membrane analysis have the same distributions. Where the macros could be improved, the
port does not reproduce them exactly. The differences are listed below, and described in detail in `docs/SPEC.md` in
the repository.

## Differences

| Topic | The macros | This program | Why |
|---|---|---|---|
| Input | a Combiner stack | the IDEAS zip, with the object numbers kept | each cell can be traced back to IDEAS, and there is one step less |
| Cellpose montage | one row of all events | a grid that is about square | faster, with the same masks |
| Cells that cross the border between two events in the montage | one cell | split into one cell per event | each event is a separate image |
| Names | T cell | immune cell | the immune side can be a panel of markers |
| Immune markers | one (1.8); raw values summed (2.8) | several, each normalized, then summed | markers of different brightness get equal weight |
| Clustering | 1.8: normalized markers and area; 2.8: raw sums, with the functional markers added without normalization | as in 1.8 | 2.8 added the functional markers without normalizing them, which was a mistake |
| Subtype thresholds | MaxEntropy on all cells (2.8) | *Two populations* (or MaxEntropy or Otsu) on the immune cells only | tumor cells no longer shift the threshold, and MaxEntropy called noise positive or missed positive cells (see [Cell types](../analysis/classification.md#immune-subtypes)) |
| Detections that are not cells | a "crap" class, switched off and not working | artefacts: detections that are positive in all subtype markers | — |
| Contacts | region adjacency graph of the labels | cells within the interface width of each other | consistent with the interfaces |
| Contacts compared with chance | — | partner preference in events with 3 or more cells, and subset preference among immune cells, with intervals | new |
| Membranes | one ring per blob of touching cells of the same class | one membrane per cell | the cell is the unit of analysis |
| Straightening width | 5 px | about 4.5 µm (9 px at 40×) | 5 px at 40× gave a false enrichment at every contact (see [Membrane enrichment](../analysis/membrane.md#the-membrane-profile)) |
| Baseline of the enrichment | all of the membrane outside the interface | the membrane that touches no other cell | in the macros, contacts with other marked cells raise the baseline |
| Control | one window of 25 px per membrane, unseeded | 10 windows as long as the median interface, on free membrane | a better null distribution |
| Average profile | all rings, including those without an interface | only cells with an interface | rings without an interface halve the peak |
| Pair of cell types | the cell of interest against the other class | any pair, including subtypes, with a check of where the marker is expressed | — |
| Several images | concatenated (2.8), with the background of the last image | pooled per experiment, with each dataset's own background | 2.8 used the background of the last image for all images |
| Record of the settings | `settings.txt`, incomplete | the complete `settings.toml` | — |

??? note "Technical: what is kept exactly"

    Several steps are the same as in macro 1.8, so that the results stay comparable: the normalization of the markers
    (the median over the events of the Otsu threshold of each event, computed with ImageJ's 8-bit histogram of float
    images), the brightfield variance filter (with ImageJ's circular kernel), the erosion and dilation (CLIJ2's
    4-connected steps), the k-means features and the rule that picks the tumor cluster, the per-event table, and the
    doublet table. The ImageJ automatic threshold methods are ported exactly and give the same thresholds as Fiji.
