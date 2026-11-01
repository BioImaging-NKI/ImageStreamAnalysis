# The project overview

The last step of the import page shows one row per analysed dataset, below it the partner preference of every
experiment and, when the datasets have a condition and a replicate, the comparison of the conditions. The same is
written to `overview.tsv`, `interactions.tsv` and `condition comparison.tsv` in the project folder.

![The project overview](../images/gui_overview.png)

## The table

One row per dataset; *dashboard* opens the dataset's dashboard, *folder* its folder. A dataset that failed shows its
error instead. Hovering a column header shows what it means:

| Column | Meaning |
|---|---|
| *events*, *cells* | the events analysed; the cells after the clean-up, artefacts included |
| *immune cells*, *tumor cells*, *artefacts* | the cell types (see [Cell types](../analysis/classification.md)) |
| *single tumor + single immune cell* | events with exactly one tumor and one immune cell that touch |
| *interactions* | pairs of touching cells, of any types |
| *MOI partners / measured* | the marker check (see [Membrane enrichment](../analysis/membrane.md#where-the-marker-is-expressed)); empty without a marker of interest |
| *marker check* | the marker check of the project's membrane pair, in words |
| *tumor–immune contacts* | pairs of a tumor and an immune cell that touch |
| *doublet fraction* | *single tumor + single immune cell* divided by *events* |
| *median Tumor (or Immune) membrane ratio (doublets)* | in those events: the mean marker of interest in the measured cell's membrane band at the interface, over the rest of its band (the macro's ratio; all immune cells, whatever the membrane pair; no control, see [the per-event table](../analysis/membrane.md#the-per-event-table)) |
| *immune negative*, *immune A+*, *immune A+ B+*, … | immune cells with exactly that combination of positive subtype markers (*immune A+* is negative for the other markers); 0 where a dataset has none |
| *membrane pair* | the pair of the enrichment columns that follow: *Membrane of* … *with* … |
| *membranes*, *with interface* | the measured cells of that pair, and those that touch a partner |
| *median stdDevs above mean*, *median control (windows)* | the median enrichment of the cells with an interface, and of the control windows |

The enrichment columns follow the membrane pair of the project: after *Apply to the project…* with another pair in a
dashboard, they describe that pair. Note the pair with the numbers. Every other pair of every dataset is in
`membrane pairs.tsv` in the project folder, without another run.

## Partner preference per experiment

One matrix per experiment, for all its datasets together: contacts observed ÷ expected by chance, in the events of
three or more cells (blue: fewer than chance, red: more). Hover a square for the counts. Comparing the matrices of
experiments shows whether a preference holds across samples. Between immune subtypes, read the
[subset preference](../analysis/interactions.md#subset-preference-which-immune-cells-touch) instead: this matrix also
compares cells of different size (immune–immune ratios well below 1 appear without any preference), and it pools the
datasets of the experiment, so it cannot show a difference between conditions within it. How it is computed, and what it does not measure, is in
[Cell contacts](../analysis/interactions.md).

## Comparing conditions

When the datasets have a *Condition* and a *Replicate* (step 2), this section tests, across the replicates, whether
the subset preference and the membrane enrichment differ between conditions: *Compare with* chooses the reference
condition, the notes say how the conditions were compared (paired or unpaired) and which datasets were left out, and a
table per measure gives the change with its 95 % interval, *p* and *q* (corrected for the number of pairs), and the
value of every replicate. Read *q*, and look whether the change goes the same way in every replicate. How to label the
datasets, how the test works and how much it can find with a few donors: [Comparing
conditions](../analysis/comparing-conditions.md).

??? note "Technical: summing over datasets"

    `interactions.tsv` holds every dataset's contact table, and per experiment a row "(all)": the contacts, the cells
    touching and the expected contacts summed over its datasets, and observed / expected computed from those sums
    (not the mean of the datasets' ratios, which would weigh a small dataset as much as a large one).
