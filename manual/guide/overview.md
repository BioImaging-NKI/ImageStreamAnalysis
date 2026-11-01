# The project overview

The last step of the import page shows one row per analysed dataset, and below it the partner preference of every
experiment. The same is written to `overview.tsv` and `interactions.tsv` in the project folder.

![The project overview](../images/gui_overview.png)

## The table

Per dataset: the events and cells; immune cells, tumor cells and artefacts; events with one tumor and one immune cell;
the pairs of touching cells and the tumor–immune contacts; the marker-of-interest check (*MOI partners / measured*,
see [Membrane enrichment](../analysis/membrane.md#where-the-marker-is-expressed)); the number of immune cells per
subtype; and with the membrane analysis the number of membranes, those with an interface, and the medians of the
enrichment and of the control. *dashboard* opens the dataset's dashboard, *folder* its folder. A dataset that failed
shows its error instead.

## Partner preference per experiment

One matrix per experiment, for all its datasets together: contacts observed ÷ expected by chance, in the events of
three or more cells (blue: fewer than chance, red: more). Hover a square for the counts. Comparing the matrices of
experiments shows whether a preference holds across samples. How it is computed, and what it does not measure, is in
[Cell contacts](../analysis/interactions.md).

??? note "Technical: summing over datasets"

    `interactions.tsv` holds every dataset's contact table, and per experiment a row "(all)": the contacts, the cells
    touching and the expected contacts summed over its datasets, and observed / expected computed from those sums
    (not the mean of the datasets' ratios, which would weigh a small dataset as much as a large one).
