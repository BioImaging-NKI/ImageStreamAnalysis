# Comparing conditions

Does a treatment change which immune cells touch the tumor cells, or how strongly a marker gathers at the
interface? The overview answers this with a test across **replicates**: the independent units of the experiment,
usually donors. It compares two measures between conditions:

- the [subset preference](interactions.md#subset-preference-which-immune-cells-touch): observed / expected contacts of
  each immune subtype with tumor cells and with the other subtypes;
- the [membrane enrichment](membrane.md): the median enrichment (SD above mean) of every membrane pair of the cell of
  interest with the other cell type.

## Labelling the datasets

In the dataset tables of the import page (step 2), give every dataset that takes part:

- a **Condition**: *control*, *treated*, … Datasets without a condition take no part.
- a **Replicate**: the donor, or the independent experiment, that the dataset comes from. Give the control and the
  treated sample of one donor the **same** replicate: the comparison is then made within each donor (paired).
  Datasets without a replicate are left out (the overview says which).

Several datasets with the same condition and replicate (for example two tubes of one donor) are pooled into one.
They are technical repeats, not replicates: only independent units count as replicates.

The labels can be set or changed after the analysis: the overview follows them at once, and *Save* writes
`condition comparison.tsv` again. From the command line, `imagestream-cli overview <project folder>` does the same.

**Classify the conditions together.** A subtype can differ between conditions by its threshold alone. Give the
datasets that are compared one experiment label (so they are [classified together](classification.md#several-datasets-pools)),
or set the thresholds by hand. The overview warns when the compared datasets were classified apart.

## Reading the result

The overview's section *Comparing conditions* compares every condition with the one chosen under *Compare with* (the
reference; by default the first condition in the project). It first says how it compared:

- **paired**, when two or more replicates have both conditions (replicates with only one of them are left out);
- **unpaired**, when the conditions come from different replicates, two or more each;
- **no test**, otherwise: the values per replicate are shown, without a test.

Then one table per measure, with the pairs that could be tested, most significant first:

| Column | Meaning |
|---|---|
| *ratio (condition / reference)* | subset preference: observed / expected in the condition divided by that in the reference; above 1: these cells touch more often after the treatment |
| *difference (median SD)* | membrane enrichment: the median enrichment in the condition minus that in the reference |
| *95% interval* | where the true change lies, with 95 % confidence, judged from the spread between the replicates |
| *p* | the chance of a change at least this large between the replicates if the conditions did not differ |
| *q* | *p* corrected for the number of pairs tested; a pair with *q* below 0.05 is shown in bold |
| *Per replicate* | the value of each replicate, reference → condition: look whether the change goes the same way in each |
| *Note* | why a pair was not tested, the marker check of a membrane pair, a warning about the classification |

A pair is not tested when any replicate has fewer than 5 contacts expected (subset preference) or fewer than 10 cells
with an interface (membrane enrichment): `condition comparison.tsv` lists these pairs too, with the reason.

**Read *q*, not *p*.** With ten or more pairs, one of them usually has *p* below 0.05 by chance alone. *q* holds the
chance of any false finding near 5 %.

**Few replicates, little power.** The test judges a change by how consistent it is between the replicates, so two
replicates rarely give a significant result, even for a real change (the interval is then wide: one degree of
freedom). Changes of contact preference are small, because few events give an immune cell a choice of partner: in
simulations of such a change, three donors with 300 events each found it in about one experiment in ten, and four
donors with 3000 events in three out of four (see below). A non-significant result with few replicates is no evidence
that nothing changed: look at the interval, which shows which changes are still possible.

!!! warning "Why there is no test per donor"
    Each dataset's subset preference has a 95 % interval from the sampling of its events. Two samples of one donor
    whose intervals do not overlap still do not show that the treatment acted: two samples of the same cells differ
    by more than the sampling of their events (preparation, staining, the gate in IDEAS). In simulations, a test of one
    donor's two samples against each other called a difference in up to three out of four comparisons where there
    was none. The tool therefore tests only across replicates, where that variation is part of the error.

??? note "Technical: the comparison, exactly"

    **Units.** A replicate–condition combination is one unit. For the subset preference the observed contacts *O* and
    the expected contacts *E* of a pair (subtype *S*, partner) are summed over the unit's datasets; its value is
    $y = \ln\big((O + 0.5) / (E + 0.5)\big)$, the log of observed / expected with 0.5 added to both so that a pair
    without contacts still has a value. Each immune–immune pair is tested once (*S*–*T* and *T*–*S* are the same
    contacts). For the membrane enrichment the value is the median SD above mean of the unit's datasets, weighted by
    their cells with an interface. Only the pairs of the cell of interest with the other cell type are tested
    (*Membrane of* the cell of interest or one of its subsets, *Interface with* the other type or its subsets): the
    same-type pairs are controls, not hypotheses.

    **Paired** (two or more replicates with both conditions): per replicate the difference $d_i = y_i^{\,\mathrm{cond}}
    - y_i^{\,\mathrm{ref}}$, and a one-sample t-test of the mean $\bar d$ against 0 with $n - 1$ degrees of freedom.
    The interval is $\bar d \pm t_{0.975,\,n-1}\, s_d / \sqrt n$. For the subset preference the effect and its
    interval are shown as ratios: $e^{\bar d}$, the geometric mean over the replicates of the ratio of the two
    observed / expected. **Unpaired** (otherwise, two or more replicates per condition): Welch's t-test of the means
    of $y$, with the Welch–Satterthwaite degrees of freedom. **q**: the Benjamini–Hochberg adjustment over the pairs
    tested of one measure and one condition.

    **Eligibility.** A pair is tested only when every unit of the comparison has $E \ge 5$ (subset preference) or at
    least 10 cells with an interface (membrane enrichment); otherwise its note names the units that fall short. The
    test needs at least two replicates per side; with exactly two the t distribution has one degree of freedom and
    $t_{0.975} = 12.7$.

    **What the test assumes.** The replicates are independent, and the per-replicate differences (on the log scale for
    the subset preference) are roughly normal. Datasets of one donor labelled as different replicates break the
    first assumption and make the test too liberal. The test says whether the conditions differ in the cells that
    were exported; the [selection in IDEAS](interactions.md#what-the-selection-in-ideas-means) and the size of the
    subtypes affect both conditions, and cancel only when they are the same in both.

    **Not independent.** The subset preferences of one dataset are linked: when one subtype takes more of the contacts
    with tumor cells, the others take fewer. A real change of one subtype therefore shows as a smaller, opposite
    change of the others; read the pairs together.

??? note "Technical: calibration by simulation"

    The tests were checked on simulated experiments with a known truth (`tools/stats/calibrate_comparison.py`, in the
    repository; about 4 minutes). Events were made with the geometry of the synthetic test data: one or two tumor
    cells and up to three immune cells of five subtypes, placed touching at random, except that a cell of subtype A
    was placed next to a tumor cell with probability $p_A$. Donors differed ($p_A$ = 0.4 with SD 0.1 between donors);
    150–3000 events per sample. Only the contacts were simulated, not the images. In this geometry $p_A$ = 0.4 gives
    a true A–tumor preference of 1.07, and $p_A$ = 0.95 one of 1.12.

    **No difference between the conditions** (proportion of experiments with *p* < 0.05; 2000 each, so ±0.01):

    | donors | paired | unpaired (Welch) | paired, samples differ (SD 0.05 on $y$) |
    |---|---|---|---|
    | 2 | 0.053 | 0.024 | 0.057 |
    | 3 | 0.049 | 0.030 | 0.043 |
    | 4 | 0.044 | 0.036 | 0.059 |
    | 6 | 0.049 | 0.047 | 0.044 |

    The paired test holds its 5 %, also when samples differ beyond the sampling of their events; Welch's test is
    conservative with few replicates.

    **A test of one donor's two samples** against each other ($z$ = difference / its standard error from the events,
    by the jackknife over events), no difference, proportion with $|z| > 1.96$:

    | SD between two samples of a donor, on $y$ | 300 events | 3000 events |
    |---|---|---|
    | 0 | 0.060 | 0.048 |
    | 0.03 | 0.066 | 0.364 |
    | 0.05 | 0.178 | 0.512 |
    | 0.10 | 0.398 | 0.754 |

    Correct only when two samples of a donor differ by nothing but their events, which cannot be checked; the more
    events, the worse. This is why the tool has no test per donor.

    **All subset pairs at once** (no difference anywhere), the chance of at least one finding:

    | donors, events per sample | pairs tested (mean) | any *p* < 0.05 | any *q* < 0.05 |
    |---|---|---|---|
    | 3, 300 | 5.3 | 0.20 | 0.053 |
    | 4, 1000 | 14.7 | 0.51 | 0.058 |
    | 6, 300 | 4.3 | 0.19 | 0.049 |

    **Power**: subtype A placed next to a tumor cell with probability 0.4 in the control and 0.95 after treatment
    (true preference 1.07 → 1.12), proportion of experiments with *p* < 0.05:

    | events per sample | 2 donors | 3 | 4 | 6 | 10 |
    |---|---|---|---|---|---|
    | 300 | 0.06 | 0.11 | 0.16 | 0.28 | 0.50 |
    | 1000 | 0.09 | 0.23 | 0.39 | 0.65 | 0.93 |
    | 3000 | 0.18 | 0.47 | 0.74 | 0.96 | 1.00 |

    A cell type placed next to the tumor cell more than twice as often changes observed / expected by only 4 %:
    most events offer little choice. Plan the number of donors and events with this in mind.

    **Membrane enrichment**: per-cell scores drawn skewed (a normal core with a quarter of the cells in an
    exponential tail), 10–200 cells with an interface per sample, donors differing; the paired test of the medians
    gave 0.043–0.053 with 2–6 donors.
