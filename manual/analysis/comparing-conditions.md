# Comparing conditions

Does a treatment change which immune cells touch the tumor cells, or how strongly a marker gathers at the
interface? The overview answers such questions with a test across **replicates**, the independent units of the
experiment (usually donors). It compares two measures between conditions:

- the [subset preference](interactions.md#subset-preference-which-immune-cells-touch): the ratio observed / expected
  of the contacts of each immune subtype with tumor cells and with the other subtypes;
- the [membrane enrichment](membrane.md): the median enrichment (in SD above the mean) of every membrane pair of the
  cell of interest with the other cell type.

## Labeling the datasets

In the dataset tables of the import page (step 2), give each dataset that should take part in the comparison:

- a **Condition**, for example *control* or *treated*. Datasets without a condition do not take part.
- a **Replicate**: the donor, or the independent experiment, that the dataset comes from. Give the control sample and
  the treated sample of one donor the **same** replicate. The comparison is then made within each donor (paired).
  Datasets without a replicate are left out, and the overview lists them.

Several datasets with the same condition and the same replicate (for example two tubes from one donor) are pooled
into one. They are technical repeats, not replicates: only independent units count as replicates.

You can set or change the labels after the analysis. The overview updates at once, and *Save* writes
`condition comparison.tsv` again. On the command line, `imagestream-cli overview <project folder>` does the same.

**Classify the conditions together.** Otherwise a subtype can differ between conditions only because its threshold
differs. Give the datasets that you compare one experiment label, so that they are
[classified together](classification.md#several-datasets-pools), or set the thresholds by hand. The overview warns
you when the datasets you compare were classified separately.

## Reading the result

The section *Comparing conditions* of the overview compares each condition with the reference condition, which you
choose under *Compare with*. By default, the reference is the first condition in the project. The section first
states how the comparison was made:

- **paired**, when two or more replicates have both conditions. Replicates that have only one of the two conditions
  are left out.
- **unpaired**, when the conditions come from different replicates, with two or more replicates for each condition.
- **no test**, in all other cases. The values of each replicate are then shown without a test.

Below this, the section has one table per measure. Each table lists the pairs that could be tested, the most
significant first:

| Column | Meaning |
|---|---|
| *ratio (condition / reference)* | For the subset preference: observed / expected in the condition, divided by observed / expected in the reference. A value above 1 means that these cells touch more often after the treatment. |
| *difference (median SD)* | For the membrane enrichment: the median enrichment in the condition minus the median enrichment in the reference. |
| *95% interval* | The range that contains the true change with 95 % confidence, judged from the spread between the replicates. |
| *p* | The probability of a change at least this large between the replicates if the conditions did not differ. |
| *q* | *p* corrected for the number of pairs tested. A pair with *q* below 0.05 is shown in bold. |
| *Per replicate* | The value of each replicate, as reference → condition. Check whether the change goes in the same direction in every replicate. |
| *Note* | Why a pair was not tested, the result of the marker check of a membrane pair, or a warning about the classification. |

A pair is not tested when any replicate has fewer than 5 expected contacts (subset preference) or fewer than 10 cells
with an interface (membrane enrichment). `condition comparison.tsv` also lists these pairs, with the reason.

**Read *q*, not *p*.** When ten or more pairs are tested, one of them usually has a *p* below 0.05 by chance alone.
*q* keeps the probability of any false finding near 5 %.

**Few replicates give little power.** The test judges a change by how consistent it is between the replicates, so
two replicates rarely give a significant result, even when the change is real. With two replicates the interval is
wide, because the test has only one degree of freedom. Changes in contact preference are small, because few events
give an immune cell a choice of partner. In simulations of such a change, experiments with three donors and 300
events per sample detected it in about one experiment in ten, and experiments with four donors and 3000 events per
sample in three out of four (see the calibration below). With few replicates, a result that is not significant is no
evidence that nothing changed. Look at the interval instead: it shows which changes are still possible.

!!! warning "Why there is no test per donor"
    The subset preference of each dataset has a 95 % interval, which reflects the sampling of its events. Even when
    the intervals of two samples from one donor do not overlap, this does not show that the treatment had an effect.
    Two samples of the same cells differ by more than the sampling of their events, for example through the
    preparation, the staining or the gate in IDEAS. In simulations, a test of one donor's two samples against each
    other reported a difference in up to three out of four comparisons in which there was none. The tool therefore
    tests only across replicates, where this variation is part of the error.

??? note "Technical: the comparison, exactly"

    **Units.** Each combination of a replicate and a condition is one unit. For the subset preference, the observed
    contacts *O* and the expected contacts *E* of a pair (subtype *S* and a partner) are summed over the datasets of
    the unit. The value of the unit is $y = \ln\big((O + 0.5) / (E + 0.5)\big)$: the log of observed / expected,
    with 0.5 added to both so that a pair without contacts still has a value. Each immune–immune pair is tested
    once, because *S*–*T* and *T*–*S* are the same contacts. For the membrane enrichment, the value of the unit is
    the median SD above mean of its datasets, weighted by their numbers of cells with an interface. Only the pairs of
    the cell of interest with the other cell type are tested (*Membrane of* the cell of interest or one of its
    subsets, and *Interface with* the other cell type or one of its subsets). The pairs within one cell type are
    controls, not hypotheses.

    **Paired** (two or more replicates with both conditions): for each replicate the difference is $d_i =
    y_i^{\,\mathrm{cond}} - y_i^{\,\mathrm{ref}}$, and a one-sample t-test compares the mean $\bar d$ with 0, with
    $n - 1$ degrees of freedom. The interval is $\bar d \pm t_{0.975,\,n-1}\, s_d / \sqrt n$. For the subset
    preference, the effect and its interval are shown as ratios: $e^{\bar d}$ is the geometric mean, over the
    replicates, of the ratio of the two observed / expected values. **Unpaired** (in all other cases, with two or more
    replicates per condition): Welch's t-test compares the means of $y$, with the Welch–Satterthwaite degrees of
    freedom. **q**: the Benjamini–Hochberg adjustment over the pairs tested for one measure and one condition.

    **Eligibility.** A pair is tested only when every unit of the comparison has $E \ge 5$ (subset preference) or at
    least 10 cells with an interface (membrane enrichment). Otherwise, its note names the units that fall short. The
    test needs at least two replicates on each side. With exactly two, the t distribution has one degree of freedom,
    and $t_{0.975} = 12.7$.

    **What the test assumes.** The replicates are independent, and the differences per replicate (on the log scale
    for the subset preference) are roughly normally distributed. Datasets from one donor that are labeled as
    different replicates violate the first assumption and make the test too liberal. The test tells whether the
    conditions differ in the cells that were exported. The [selection in IDEAS](interactions.md#what-the-selection-in-ideas-means)
    and the size of the subtypes affect both conditions, and their effects cancel only when they are the same in
    both.

    **Not independent.** The subset preferences within one dataset are linked: when one subtype takes more of the
    contacts with tumor cells, the other subtypes take fewer. A real change in one subtype therefore also shows up
    as a smaller, opposite change in the others, so read the pairs together.

??? note "Technical: calibration by simulation"

    The tests were checked on simulated experiments with a known truth (`tools/stats/calibrate_comparison.py` in the
    repository, which takes about 4 minutes). The events were made with the geometry of the synthetic test data:
    one or two tumor cells and up to three immune cells of five subtypes, placed touching each other at random,
    except that a cell of subtype A was placed next to a tumor cell with probability $p_A$. The donors differed
    ($p_A$ = 0.4, with an SD of 0.1 between donors), and each sample had 150–3000 events. Only the contacts were
    simulated, not the images. In this geometry, $p_A$ = 0.4 gives a true A–tumor preference of 1.07, and
    $p_A$ = 0.95 gives one of 1.12.

    **No difference between the conditions** (the proportion of experiments with *p* < 0.05; 2000 experiments each,
    so the proportions are accurate to ±0.01):

    | donors | paired | unpaired (Welch) | paired, samples differ (SD 0.05 on $y$) |
    |---|---|---|---|
    | 2 | 0.053 | 0.024 | 0.057 |
    | 3 | 0.049 | 0.030 | 0.043 |
    | 4 | 0.044 | 0.036 | 0.059 |
    | 6 | 0.049 | 0.047 | 0.044 |

    The paired test keeps its error rate at 5 %, also when the samples differ by more than the sampling of their
    events. Welch's test is conservative when there are few replicates.

    **A test of one donor's two samples** against each other ($z$ = the difference divided by its standard error
    from the events, estimated by the jackknife over events), with no true difference; the table shows the
    proportion of comparisons with $|z| > 1.96$:

    | SD between two samples of a donor, on $y$ | 300 events | 3000 events |
    |---|---|---|
    | 0 | 0.060 | 0.048 |
    | 0.03 | 0.066 | 0.364 |
    | 0.05 | 0.178 | 0.512 |
    | 0.10 | 0.398 | 0.754 |

    This test is correct only when two samples of a donor differ in nothing but their events, which cannot be
    checked. The more events, the worse it gets. This is why the tool has no test per donor.

    **All subset pairs at once** (with no difference anywhere), the probability of at least one finding:

    | donors, events per sample | pairs tested (mean) | any *p* < 0.05 | any *q* < 0.05 |
    |---|---|---|---|
    | 3, 300 | 5.3 | 0.20 | 0.053 |
    | 4, 1000 | 14.7 | 0.51 | 0.058 |
    | 6, 300 | 4.3 | 0.19 | 0.049 |

    **Power**: subtype A was placed next to a tumor cell with probability 0.4 in the control and 0.95 after the
    treatment (a true preference of 1.07 → 1.12). The table shows the proportion of experiments with *p* < 0.05:

    | events per sample | 2 donors | 3 | 4 | 6 | 10 |
    |---|---|---|---|---|---|
    | 300 | 0.06 | 0.11 | 0.16 | 0.28 | 0.50 |
    | 1000 | 0.09 | 0.23 | 0.39 | 0.65 | 0.93 |
    | 3000 | 0.18 | 0.47 | 0.74 | 0.96 | 1.00 |

    A cell type that is placed next to the tumor cell more than twice as often changes observed / expected by only
    4 %, because most events offer little choice. Keep this in mind when you plan the number of donors and events.

    **Membrane enrichment**: the scores per cell were drawn from a skewed distribution (a normal core, with a quarter
    of the cells in an exponential tail), with 10–200 cells with an interface per sample and donors that differed.
    With no difference between the conditions, the paired test of the medians gave *p* < 0.05 in 0.043–0.053 of the
    experiments, with 2–6 donors.
