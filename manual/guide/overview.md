# The project overview

The last step of the import page shows the project overview. It has one row per analyzed dataset. Below the table, it
shows the partner preference of every experiment and, when the datasets have a condition and a replicate, the
comparison of the conditions. The same information is written to `overview.tsv`, `interactions.tsv` and
`condition comparison.tsv` in the project folder. Only the datasets of the project count (see
[During and after a run](#during-and-after-a-run)).

![The project overview](../images/gui_overview.png)

## The table

The table has one row per dataset. The *Dashboard* button opens the dashboard of the dataset, *Folder* opens its
folder, and *Run again* runs only that dataset again (see [Run](import-page.md#4-run)). For a dataset that failed, the row shows the error instead. Hover over a column header to see what the column
means:

| Column | Meaning |
|---|---|
| *events*, *cells* | the number of events analyzed; the number of cells after the clean-up, including artefacts |
| *immune cells*, *tumor cells*, *artefacts* | the number of cells of each type (see [Cell types](../analysis/classification.md)) |
| *single tumor + single immune cell* | the number of events with exactly one tumor cell and one immune cell that touch |
| *interactions* | the number of pairs of touching cells, of any cell types |
| *tumor–immune contacts* | the number of pairs of a tumor cell and an immune cell that touch |
| *doublet fraction* | *single tumor + single immune cell* divided by *events* |
| *immune negative*, *immune A+*, *immune A+ B+*, … | the number of immune cells with exactly that combination of positive subtype markers (*immune A+* cells are negative for the other markers); 0 when a dataset has none |

Then come the columns of each marker of interest. Their names end in the name of the marker, for example
*MOI on (Marker)*:

| Column | Meaning |
|---|---|
| *MOI on* | the cell type (*Tumor* or *Immune*) that carries the marker, found in the data; by default its membranes are measured (see [A pair of cell types](../analysis/membrane.md#a-pair-of-cell-types)) |
| *MOI partners / measured* | the marker check (see [Membrane enrichment](../analysis/membrane.md#where-the-marker-is-expressed)) |
| *marker check* | the marker check of the marker's membrane pair, in words |
| *median Tumor (or Immune) membrane ratio, doublets* | computed in the events of *single tumor + single immune cell*, for the cell type that carries the marker: the mean of the marker in that cell's membrane band at the interface, divided by that in the rest of its band. This is the macro's ratio. It uses all partners of the other type, whatever the membrane pair, and has no control; see [the per-event table](../analysis/membrane.md#the-per-event-table) |
| *membrane pair* | the pair that the enrichment columns after it describe: *Membrane of* … *with* … |
| *membranes*, *with interface* | the number of measured cells of that pair, and how many of them touch a partner |
| *median log2 ratio*, *median control log2 ratio (windows)* | the median enrichment of the cells with an interface, as the log2 ratio of the interface and the free membrane (1 = twice as bright), and that of the control windows (about 0); the main readout |
| *median stdDevs above mean*, *median control (windows)* | the same in standard deviations of the free membrane (the macro's score) |

The enrichment columns describe the membrane pair of each marker. If you click *Apply to the project…* in a dashboard
with another pair selected, they describe that pair, so always note the pair together with the numbers. Every other
pair of every dataset is in `membrane pairs.tsv` in the project folder, so you don't need another run to see them.

## During and after a run

The overview reads the results that are in the project folder. While a run is going on, some of them are about to be
replaced, so a column before *experiment* shows the state of each row:

| State | Meaning |
|---|---|
| *running…* (greyed out) | the run is computing this dataset again; until it is done, the values are those of the previous run |
| *running…* (greyed out, no links) | the run is computing this dataset for the first time; there are no results yet |
| *updated* | written by the last run |
| *previous run* (greyed out) | the last run was stopped or failed before it reached this dataset; the values come from an earlier run, with the settings of that time |
| *not in the project* (greyed out) | see below |

The overview is refreshed every time a dataset is done, and once more at the end. The partner preference and the
comparison of conditions are computed at the end of a run, so during a run they are greyed out with a note.

**Results that are not part of the project.** If you reuse an output folder, it can hold result folders of datasets
that are no longer in the project: a dataset that was removed or unticked in step 2, or one whose experiment label was
changed (its results go to a folder named after the experiment). These results are left out of the overview, of
`overview.tsv`, of the partner preference and of the comparison. A line under the table says how many there are, and
*Show them* lists them, greyed out and marked *not in the project*. Their files stay on disk: delete those folders
yourself if you no longer need them.

## Partner preference per experiment

For each experiment, a matrix shows the contacts observed divided by the contacts expected by chance, for all its
datasets together. Only events with three or more cells are used. Blue means fewer contacts than chance, and red means
more. Hover over a square to see the counts. Comparing the matrices of different experiments shows whether a
preference holds across samples.

Between immune subtypes, read the
[subset preference](../analysis/interactions.md#subset-preference-which-immune-cells-touch) instead, for two reasons.
First, this matrix also compares cells of different sizes, so immune–immune ratios well below 1 can appear without any
preference. Second, it pools the datasets of the experiment, so it cannot show a difference between conditions within
that experiment. [Cell contacts](../analysis/interactions.md) explains how the matrix is computed and what it does not
measure.

## Comparing conditions

When the datasets have a *Condition* and a *Replicate* (set in step 2), this section tests whether the subset
preference and the membrane enrichment differ between conditions, across the replicates. Use *Compare with* to choose
the reference condition. The notes say how the conditions were compared (paired or unpaired) and which datasets were
left out. For each measure, a table gives the change with its 95 % interval, *p* and *q* (corrected for the number of
pairs tested), and the value of every replicate. Read *q*, and check whether the change goes in the same direction in
every replicate. [Comparing conditions](../analysis/comparing-conditions.md) explains how to label the datasets, how
the test works, and how much it can detect with a few donors.

??? note "Technical: summing over datasets"

    `interactions.tsv` holds the contact table of every dataset and, for each experiment, a row "(all)". In that row,
    the contacts, the touching cells and the expected contacts are summed over the datasets of the experiment, and
    observed / expected is computed from those sums. It is not the mean of the ratios of the datasets, because that
    would give a small dataset as much weight as a large one.
