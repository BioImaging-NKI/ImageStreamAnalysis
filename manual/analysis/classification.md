# Cell types

Every detection is classified as a **tumor cell**, as an **immune cell** with a **subtype**, or as an **artefact**.
This happens in two steps. The first step decides the lineage: tumor or immune. The second step assigns a subtype to
each immune cell, based on the subtype markers.

## Tumor or immune cell

The cells are split into two groups by k-means clustering on three values per cell: the tumor marker and the immune
marker (both normalized and on a log scale), and the area. The group with the lower ratio of immune to tumor marker is
taken to be the tumor cells. This is the classification of the Fiji macro v1.8.

??? note "Technical: the clustering"

    The features are log(immune intensity), log(tumor intensity) and area, each scaled to a mean of 0 and an SD of 1
    (using the sample SD). The k-means clustering uses 2 clusters, 10 starts and seed 48, which is the seed of the
    macro. The k-means of scikit-learn is not bit-identical to that of ImageJ, but on well-separated data the result
    is the same: the classes agreed for 100 % of the cells on the validation data. The tumor cluster is the one with
    the lower **sum** of *marker ratio over area*, as the macro computes it. (The comment in the macro says
    "average"; with clusters of unequal size, the sum and the average can give different results.)

    When there are several immune markers, they are summed after normalization. The immune cells can therefore be
    defined by a panel of markers (for example CD3 and CD56) instead of by a single marker.

    A cell whose mean tumor or immune intensity is 0 or less has no logarithm. This can happen when part of a channel
    of an event image is empty (all pixels 0) where the cell lies, a defect of the exported image. Such a cell gets
    half the lowest positive intensity of the other cells instead, so that it is clustered by its other marker, and
    the log names its event.

## Immune subtypes

Each **subtype marker** gets a threshold, and an immune cell above the threshold is positive for that marker. The
combination of markers for which a cell is positive is its subtype, and the subtype is named after those markers: for
example *Subtype A+*, *Subtype A+ Subtype B+*, or *negative* when the cell is positive for none of them. The thresholds
are either found automatically or set by hand, in the settings (`classification.thresholds`) or with the sliders in
the [dashboard](../guide/dashboard.md#subtype-markers).

The setting *Threshold of the subtype markers* chooses the automatic method:

- **Two populations** (the default): looks for a second, positive population next to the negative cells, and puts the
  threshold where a cell is equally likely to belong to either population. If there is no such population (because
  the positive cells form a tail rather than a peak), the threshold is the upper noise limit of the negative cells.
  This method handles a small positive population close to a large negative one, which is the usual situation for a
  subtype marker.
- **MaxEntropy** (as in ClassifyCellTypes 2.8): finds the split of the histogram that gives the most information on
  both sides.
- **MaxEntropy strict**: the MaxEntropy threshold, raised where necessary so that it lies at least two robust standard
  deviations above the middle of the negative cells. This gives fewer false positives, at the cost of a few more
  missed positives.
- **Otsu**: the split that best separates two populations of similar size.

For every threshold, the `how` column of `Marker thresholds.tsv` states how it was found.

**Check every threshold on its histogram** in the dashboard before you interpret the subtype numbers. Where the
negative and positive cells form two clear peaks, the threshold belongs in the gap between them. The automatic methods
can miss this gap on either side (see the technical note). If you set the threshold in the gap with the slider and
apply it to the project (*Apply to the project…*), it is then used for every output.

![The subtype markers in the dashboard](../images/dash_markers.png){ width="420" }

??? note "Technical: the thresholds"

    Marker intensities are roughly log-normally distributed, so the threshold is found on

    $$x = \ln(v - \min v + 1)$$

    where *v* is the background-subtracted mean intensity of a cell. The threshold $t_x$ is found with ImageJ's
    auto-threshold on a 256-bin histogram of *x* (as `setAutoThreshold("MaxEntropy dark")` in 2.8), and then
    transformed back: $t = e^{t_x} - 1 + \min v$. A cell is positive when *v* > *t*. The MaxEntropy, Otsu, Yen and
    RenyiEntropy methods are ported from ImageJ's AutoThresholder and give the same thresholds as Fiji (checked to
    10⁻¹¹).

    **The thresholds are computed on the immune cells only.** A subtype marker may also be expressed by tumor cells.
    Computing its threshold on the immune cells alone prevents the tumor cells from shifting it.

    **The automatic thresholds can be wrong in either direction.** MaxEntropy often cuts in the upper tail of the
    negative population, so that some negative cells are called positive. On synthetic data this typically happens to
    10–20 % of them, or about ten cells per marker in a dataset of a few hundred immune cells. This creates small,
    spurious subtypes such as *A+ C+*. On other data, MaxEntropy cuts into the positive population instead: in a test
    at 60×, it missed a quarter of the positive cells, all of which Otsu found. *MaxEntropy strict* can only raise the
    threshold. The thresholds are listed in `Marker thresholds.tsv`, as background-subtracted intensities (the same
    scale as `Mean ch<N>`).

??? note "Technical: Two populations, step by step"

    All steps work on the log scale *x* defined above, with one value per cell of the population: the immune cells,
    or, for the artefact rule, the immune cells of the first clustering.

    1. **The negative cells.** A kernel density estimate of *x* is computed (Gaussian kernel; bandwidth
       $0.9 \cdot \min(s, \mathrm{IQR}/1.34) \cdot n^{-1/5}$ by Silverman's rule; evaluated at 1024 points). Its
       highest peak is taken as the mode *m* of the negative cells. Their spread is estimated from the **lower** side
       of the peak only, which positive cells do not reach: $\sigma = 1.4826 \cdot \mathrm{median}(m - x)$ over the
       cells with $x \le m$. For a normal distribution, this equals the standard deviation. The **noise limit** is
       $m + 2.5\,\sigma$.
    2. **A second population.** A mixture of two normal distributions is fitted to *x* by expectation–maximization,
       starting from the mode and from the largest value, so that the result is deterministic. The upper component
       counts as a positive population when it contains at least 0.5 % of the cells, when its mean lies at least two
       pooled standard deviations from the mean of the lower component (Ashman's *D* ≥ 2), and when its mean lies
       above the noise limit.
    3. **The threshold.** If there is a positive population, the threshold is where the two weighted normal densities
       cross, between their means: above this point, a cell is more likely to be positive than negative. The threshold
       is never lower than $m + 2\sigma$. If there is no positive population, the threshold is the noise limit. With
       fewer than ten cells, MaxEntropy is used instead.

    The `how` column of `Marker thresholds.tsv` reports which case applied, together with the fraction of cells in
    the fitted positive population: *two populations (12 % positive): crossing*, *... raised to mode + 2 SD*, or
    *no second population: noise limit*.

    **Why this method.** A threshold that maximizes a criterion over the whole histogram (Otsu, MaxEntropy) is pulled
    into the large negative population when the positive population is small, which is the common situation for a
    subtype marker. A fixed distance above the negative cells (as in "strict") is right when there is no clear second
    population, but it cuts into a positive population that lies close to the negative one. The two steps cover both
    situations. The method was chosen on simulated data: 600 cells, with log-normal negative cells with noise (two
    spreads), positive cells making up 2–40 % of the cells and 3–40 times as bright, and five repeats (500 samples in
    all). The total number of errors (false plus missed positives), compared with the best possible threshold for each
    sample, was: best possible 2417, **Two populations 4790**, the noise limit alone 6825, MaxEntropy strict 14 881,
    Otsu 21 132 and MaxEntropy 21 965. Almost all remaining errors occur in samples where the two populations overlap
    (positive cells 2–5 % of the cells and at most 5 times as bright), where no threshold can separate them. On the
    four user-test datasets (2289 positive cells in all), Two populations made 4 errors and MaxEntropy 286.

    **Limits.** When the positive cells form a continuous tail, the noise limit decides, and the number of positive
    cells depends on the choice of 2.5 SD. Set the threshold by hand if the histogram shows a better place for it. A
    bimodal *negative* population (for example, two batches of differently stained cells pooled together) can be
    mistaken for two populations. Classify such datasets separately (*Classify together* = dataset).

??? note "Technical: why two steps, and not one clustering of all cells"

    A single k-means clustering of all cells into k types was considered and not used, for three reasons. k-means
    favors groups of similar size and spread, so a small subtype is absorbed by a large neighboring group. The number
    k must be chosen in advance. And the clusters have to be named afterwards. A threshold per marker, in contrast,
    can be read, checked on its histogram and changed, and the combinations of positive markers name the subtypes
    directly.

## Artefacts

A detection that is positive for **all** subtype markers (when there are two or more of them) is not taken to be a
cell, but an **artefact**: for example debris, or a dead cell that takes up every dye. Artefacts are left out of the
tumor/immune clustering and drawn in gray. They are not counted as cells, neither in the cell counts nor in the
interfaces. Use *Artefact: positive in at least N subtype markers* to set a different number; 0 switches the rule
off.

**The images cannot distinguish an artefact from a real cell that is positive for every marker.** If you expect such
cells, look at the artefacts in the dashboard (they are drawn in gray) before you switch the rule off. With the rule
off, real artefacts appear as a subtype that is positive for every marker; they are counted, and they take part in
the contacts as if they were cells. With **one** subtype marker, the rule is off and an artefact counts as a positive
immune cell (the log says so). In that case, check the positive cells by eye.

??? note "Technical: the artefact rule"

    The rule needs the subtype thresholds of the immune cells. These are only known after the tumor/immune
    clustering, but that clustering should run without the artefacts. The clustering therefore runs twice:

    1. A first k-means clustering of **all** detections finds the immune cells.
    2. Every subtype marker is thresholded on those immune cells, with the threshold method of the project or with the
       threshold set by hand. A detection that is positive for at least *N* markers by these thresholds is an
       artefact. The column `Positive markers` in the cell table holds the number of positive markers, and the rows
       "artefact rule" in `Marker thresholds.tsv` hold these thresholds.
    3. The final k-means clustering runs on the remaining detections. The subtype thresholds are then computed again
       on its immune cells (the rows "immune cells"), and may differ slightly.

    Until November 2026, the thresholds of the rule were computed over all detections. For well-separated markers,
    these thresholds lay at the noise level of the many negative cells (including the tumor cells). A cell positive
    for one marker that touched a cell positive for another marker picked up a little of its neighbor's light, and
    so became an artefact. In a test with three subtype markers, 5 of 47 "artefacts" were such cells. With the
    thresholds of the immune cells, the rule found exactly the 42 artefacts that had been simulated.

    Only the subtype markers count. Counting the tumor and immune markers as well (a detection that is bright for both
    could be two cells drawn as one) was tried and dropped. An automatic threshold of those markers over all
    detections calls far more cells positive than the clustering finds, so most detections would become artefacts.

## Several datasets: pools

With **Classify together** = *experiment* (the default), the cells of all datasets with the same experiment label are
classified together, with one clustering and one set of thresholds. A cell that is "Subtype A+" is then positive by
the same standard in every dataset of the experiment, so the subtype fractions can be compared between datasets. With
*dataset*, each dataset gets its own thresholds.

??? note "Technical: when to pool"

    Pool datasets that were stained and acquired in the same way (an experiment). Each dataset keeps its own
    background and its own normalization of the tumor and immune markers; only the clustering and the thresholds are
    shared. Datasets from different staining batches or with different instrument settings can differ in brightness.
    Give them different experiment labels, or classify each dataset separately. With pooling, the run first measures
    all datasets, then classifies each experiment, and then finishes each dataset. `run.json` records what each
    dataset was classified with.

## Colors

The dashboard gives up to eight cell types their own color (tumor, up to six immune types, and artefacts), in a fixed
order, so that a cell type keeps its color when other types appear or disappear. Any further immune types share a
gray color. The legends and tables always name the cell types.
