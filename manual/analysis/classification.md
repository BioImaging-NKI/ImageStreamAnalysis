# Cell types

Every detection becomes a **tumor cell**, an **immune cell** with its **subtype**, or an **artefact**, in two steps:
first the lineage, tumor or immune; then, for the immune cells, the subtype markers.

## Tumor or immune cell

The cells are split into two groups by k-means clustering on three values per cell: the tumor marker, the immune
marker (both normalized, on a log scale) and the area. The group with the lower immune/tumor marker ratio is the tumor
cells. This is the classification of the Fiji macro v1.8.

??? note "Technical: the clustering"

    Features: log(immune intensity), log(tumor intensity) and area, each scaled to mean 0 and SD 1 (sample SD).
    k-means with 2 clusters, 10 starts, seed 48 (the macro's seed; scikit-learn's k-means is not bit-identical to
    ImageJ's, but on well-separated data the result is the same: 100 % agreement on the validation data). The tumor
    cluster is the one with the lower **sum** of *marker ratio over area*, as the macro computes it (its comment says
    average; with clusters of unequal size the two can differ).

    Several immune markers are summed after normalization, so the immune side can be defined by a panel (for example
    CD3 and CD56) instead of one marker.

## Immune subtypes

Each **subtype marker** gets a threshold: an immune cell above it is positive for that marker. The combination of
positive markers is its subtype, named after them: *Subtype A+*, *Subtype A+ Subtype B+*, or *negative*. The
thresholds are automatic, or set by hand (`classification.thresholds`, or with the sliders in the
[dashboard](../guide/dashboard.md#subtype-markers)).

*Threshold of the subtype markers*:

- **Two populations** (default): looks for a second, positive population next to the negative cells and puts the
  threshold where a cell is as likely to belong to one as to the other. If there is no such population (the positive
  cells form a tail rather than a peak), the threshold is the upper noise limit of the negative cells. It copes with a
  small positive population close to a large negative one, the usual case for a subtype marker.
- **MaxEntropy** (as ClassifyCellTypes 2.8): finds the split of the histogram with the most information on both
  sides.
- **MaxEntropy strict**: the MaxEntropy threshold, raised where needed so that it lies at least two robust standard
  deviations above the middle of the negative cells. Fewer false positives, a few more missed positives.
- **Otsu**: the split that best separates two populations of similar size.

`Marker thresholds.tsv` says for every threshold how it was found (column `how`).

**Check every threshold on its histogram** in the dashboard before reading subtype numbers. Where the negative and
positive cells form two clear peaks, the threshold belongs in the gap between them; the automatic methods can miss it
on either side (see the technical note). A threshold set in the gap with the slider and applied to the project
(*Apply to the project…*) then holds for every output.

![The subtype markers in the dashboard](../images/dash_markers.png){ width="420" }

??? note "Technical: the thresholds"

    Intensities of markers are roughly log-normal, so the threshold is found on

    $$x = \ln(v - \min v + 1)$$

    of the background-subtracted cell means *v*, with ImageJ's auto-threshold on a 256-bin histogram of *x* (as
    `setAutoThreshold("MaxEntropy dark")` in 2.8), and transformed back: $t = e^{t_x} - 1 + \min v$. A cell is positive
    when *v* > *t*. The MaxEntropy, Otsu, Yen and RenyiEntropy methods are ported from ImageJ's AutoThresholder and
    give the same thresholds as Fiji (checked to 10⁻¹¹).

    **The thresholds are computed on the immune cells only.** A subtype marker may also be expressed by tumor cells;
    thresholding it within the immune cells keeps them from shifting the threshold.

    **The automatic thresholds can be off in both directions.** MaxEntropy often cuts in the upper tail of the
    negative population, so some negative cells are called positive (on synthetic data typically 10–20 % of them, or
    about ten cells per marker in a dataset of a few hundred immune cells); this creates small spurious subtypes such
    as *A+ C+*. On other data it cuts into the positive population instead (in a test at 60×, a quarter of the
    positive cells were missed, where Otsu found them all). *MaxEntropy strict* only raises the threshold. The
    thresholds are in `Marker thresholds.tsv`, in background-subtracted intensity (as `Mean ch<N>`).

??? note "Technical: Two populations, step by step"

    All on the log scale *x* above (one value per cell of the population: the immune cells; for the artefact rule the
    artefact rule).

    1. **The negative cells.** A kernel density estimate of *x* (Gaussian kernel; bandwidth
       $0.9 \cdot \min(s, \mathrm{IQR}/1.34) \cdot n^{-1/5}$, Silverman's rule; evaluated on 1024 points) has its
       highest peak at the mode *m* of the negative cells. Their spread is estimated from the **lower** side of the
       peak only, where positive cells do not reach: $\sigma = 1.4826 \cdot \mathrm{median}(m - x)$ over the cells
       with $x \le m$ (for a normal distribution this is its standard deviation). The **noise limit** is
       $m + 2.5\,\sigma$.
    2. **A second population.** A mixture of two normal distributions is fitted to *x* by expectation–maximization,
       started at the mode and at the largest value (deterministic). Its upper component counts as a positive
       population when it holds at least 0.5 % of the cells, its mean lies at least two pooled standard deviations
       from the lower one (Ashman's *D* ≥ 2), and its mean lies above the noise limit.
    3. **The threshold.** With a positive population: where the two weighted normal densities cross, between their
       means: above it a cell is more likely positive than negative. It is never lower than $m + 2\sigma$. Without
       one: the noise limit. Fewer than ten cells: MaxEntropy.

    `how` in `Marker thresholds.tsv` reports which case applied, with the fraction of the fitted positive
    population: *two populations (12 % positive): crossing*, *... raised to mode + 2 SD*, or *no second population:
    noise limit*.

    **Why.** A threshold that maximizes a criterion over the whole histogram (Otsu, MaxEntropy) is pulled into the
    large negative population when the positive one is small, which is the common case for a subtype marker; a fixed
    distance above the negative cells ("strict") is right when there is no clear second population but cuts into
    a positive population that lies close. The two steps cover both. The method was chosen on simulated data: 600
    cells, log-normal negative cells with noise (two spreads), positive cells 2–40 % of them and 3–40 times as
    bright, five repeats (500 samples). Total errors (false + missed positives) against the best possible threshold of
    each sample: best 2417, **Two populations 4790**, the noise limit alone 6825, MaxEntropy strict 14 881, Otsu
    21 132, MaxEntropy 21 965. Almost all remaining errors are in samples where the two populations overlap (positive
    cells 2–5 % and at most 5 times as bright), where no threshold separates them. On the four user-test datasets
    (2289 positive cells in all), Two populations made 4 errors, MaxEntropy 286.

    **Limits.** When the positive cells are a continuous tail, the noise limit decides, and the count depends on
    2.5 SD: set the threshold by hand if the histogram shows a clearer place. A bimodal *negative* population (two
    batches of differently stained cells pooled) can be taken for two populations: classify such datasets separately
    (*Classify together* = dataset).

??? note "Technical: why two steps, and not one clustering of all cells"

    One k-means of all cells into k types was considered and not used: k-means favours groups of similar size and
    spread, so a small subtype is absorbed by a large neighbour; k must be chosen in advance; and the clusters must be
    named afterwards. A threshold per marker can be read, checked on its histogram and changed, and the combinations
    of positive markers name the subtypes directly.

## Artefacts

A detection that is positive in **all** subtype markers (with two or more of them) is not taken to be a cell but an
**artefact**, for example debris or a dead cell that takes up every dye. Artefacts are left out of the tumor/immune
clustering, drawn in grey, and not counted as cells in the counts and the interfaces. *Artefact: positive in at least N
subtype markers* sets another number; 0 switches the rule off.

**The images cannot tell an artefact from a real cell that is positive for every marker.** If such cells are
expected, look at the artefacts in the dashboard (they are drawn in grey) before switching the rule off: with the rule
off, real artefacts appear as a subtype positive for every marker, and are counted and take part in the contacts like
cells. With **one** subtype marker the rule is off and an artefact counts as a positive immune cell (the log says
so); check the positive cells by eye.

??? note "Technical: the artefact rule"

    The rule needs the subtype thresholds of the immune cells, which are only known after the tumor/immune clustering,
    which in turn should run without the artefacts. So the clustering runs twice:

    1. a first k-means of **all** detections finds the immune cells;
    2. every subtype marker is thresholded on those immune cells (the threshold method of the project, or the
       threshold set by hand); a detection positive in at least *N* markers by these thresholds is an artefact
       (`Positive markers` in the cell table holds the count; `Marker thresholds.tsv`: the rows "artefact rule");
    3. the final k-means runs on the other detections, and the subtype thresholds are computed again on its immune
       cells (the rows "immune cells"), which may differ slightly.

    Until November 2026 the thresholds of the rule were computed over all detections. For well-separated markers they
    lay at the noise of the many negative cells (tumor cells included), so a cell positive for one marker that
    touched a cell positive for another picked up a little of its neighbour's light and became an artefact. In a test
    with three subtype markers, 5 of 47 "artefacts" were such cells; with the immune-cell thresholds the rule found
    exactly the 42 artefacts that were drawn.

    Only the subtype markers count. Counting the tumor and immune markers as well (a detection bright for both could
    be two cells drawn as one) was tried and dropped: an automatic threshold of those markers over all detections
    calls far more cells positive than the clustering finds, so most detections would become artefacts.

## Several datasets: pools

With **Classify together** = *experiment* (the default), the cells of all datasets with the same experiment label are
classified together: one clustering and one set of thresholds. A cell that is "Subtype A+" is then positive by the
same standard in every dataset of the experiment, so the subtype fractions can be compared between datasets. With
*dataset*, each dataset gets its own thresholds.

??? note "Technical: when to pool"

    Pool datasets that were stained and acquired alike (an experiment). Each dataset keeps its own background and its
    own normalization of the tumor and immune markers; only the clustering and the thresholds are shared. Datasets of
    different staining batches or instrument settings can differ in brightness: give them different experiment
    labels, or classify per dataset. The run then measures all datasets first, classifies each experiment, and
    finishes each dataset; `run.json` says what a dataset was classified with.

## Colours

The dashboard colours up to eight cell types apart (tumor, up to six immune types, artefacts), in a fixed order so that
a type keeps its colour when others come and go. More immune types share a grey; the legends and tables always name
the types.
