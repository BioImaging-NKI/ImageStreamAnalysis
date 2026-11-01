# From the Fiji macros

ImageStream Analysis for Python ports two Fiji macros of
[BioImaging-NKI/ImageStreamAnalysis](https://github.com/BioImaging-NKI/ImageStreamAnalysis):

- **`Analyze_ImageStream_stack_2celltypes_+_membrane_intensity_analysis_1.8.ijm`**: tumor and T cells by k-means, the
  per-event interface measurements, and the membrane enrichment;
- **`ClassifyCellTypes_2_8.ijm`**: several immune markers, immune subtypes by thresholds, the cells that touch, and
  a class of non-cells.

The port reads IDEAS exports directly (it replaces ImageStreamCombiner) and joins the two macros: the 1.8 analysis is
the case "one immune marker, no subtype markers". It was validated against the output of the 1.8 macro: the cells and
the per-event tables match to three decimals with 100 % class agreement, and the membrane analysis matches on its
distributions. The port does not reproduce the macros exactly where they could be improved; the differences are
below, and in detail in `docs/SPEC.md` of the repository.

## Differences

| Topic | The macros | This program | Why |
|---|---|---|---|
| Input | a Combiner stack | the IDEAS zip, object numbers kept | traceable to IDEAS; one step less |
| Cellpose montage | one row of all events | an about-square grid | faster, same masks |
| Cells across two events in the montage | one cell | split per event | each event is its own image |
| Names | T cell | immune cell | the immune side can be a panel of markers |
| Immune markers | one (1.8); summed raw values (2.8) | several, each normalized, then summed | markers of different brightness weigh alike |
| Clustering | 1.8: normalized markers and area; 2.8: raw sums, functional markers added unnormalized | 1.8's | 2.8 added the functional markers without normalizing them (a mistake) |
| Subtype thresholds | MaxEntropy on all cells (2.8) | MaxEntropy (or others) on the immune cells | tumor cells do not shift the threshold |
| Non-cells | "crap" class, switched off and not working | artefacts: positive in all subtype markers | — |
| Contacts | region adjacency graph of the labels | within the interface width | agrees with the interfaces |
| Contacts against chance | — | partner preference in events of 3+ cells | new |
| Membranes | one ring per blob of touching cells of a class | one membrane per cell | a cell is the unit |
| Baseline of the enrichment | all membrane outside the interface | the membrane that touches no cell | contacts with other marked cells raise the baseline |
| Control | one window of 25 px per membrane, unseeded | 10 windows of the median interface length, on free membrane | a better null distribution |
| Average profile | all rings, also without interface | cells with an interface | rings without interface halve the peak |
| Pair of cell types | the cell of interest against the other class | any pair, also subtypes; checked for where the marker is | — |
| Several images | concatenated (2.8), background of the last image | pools per experiment, each dataset its own background | 2.8 used the background of the last image for all |
| Settings record | `settings.txt`, incomplete | the complete `settings.toml` | — |

??? note "Technical: what is kept exactly"

    The normalization of the markers (the median over the events of each event's Otsu threshold, with ImageJ's
    8-bit histogram of float images), the brightfield variance filter (ImageJ's circular kernel), the erosion and
    dilation (CLIJ2's 4-connected steps), the k-means features and the rule that picks the tumor cluster, the
    per-event table and the doublet table are as in macro 1.8, so that results stay comparable. The ImageJ
    auto-threshold methods are ported exactly; they give Fiji's thresholds.
