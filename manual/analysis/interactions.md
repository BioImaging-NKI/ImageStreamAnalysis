# Cell contacts

## Which cells touch

Two cells of the same event **touch** when one comes within *Membrane / interface width* pixels of the other: the
same distance that defines the [interfaces](membrane.md). So a cell with an interface always lists its partner among
the cells it touches. Cells in different events never touch.

Every cell lists the cells it touches (`Touching (labels)`) and their types (`Neighbour types`); every event counts its
pairs (`Interactions`).

??? note "Technical: the distance"

    "Within *k* pixels" means *k* steps between 4-connected pixels: the reach of the dilation that makes the
    interfaces (as CLIJ2 in the macro). Two cells with no background pixel between them are 1 step apart; with a gap of
    *g* pixels, *g* + 1. ClassifyCellTypes 2.8 used MorphoLibJ's region adjacency graph on the shrunk labels, which
    counts only cells that touch directly; with the default width of 3 px, cells that almost touch also count here,
    as they do for the interfaces.

## Counting contacts

`Interactions.tsv` (and the [dashboard](../guide/dashboard.md#interactions)) count, per cell type and the type it
touches:

- **contacts**: the touching pairs; a pair of A and B counts once in row A and once in row B;
- **cells touching**: the cells of the row's type that touch at least one cell of the column's type, and that
  **fraction** of the row's cells.

## Partner preference: contacts compared with chance

Plain counts mostly reflect how many cells of each type there are. To see whether some types **prefer** each other,
the contacts are compared with chance: the contacts expected if the cell types were shuffled among the cells of each
event, with the cells of every event and its number of contacts kept as they are. The ratio **observed / expected** is
above 1 when two types touch more often than random pairing of the cells in their events would give, and below 1 when
they touch less often.

**Only events of three or more cells take part.** In an event of two touching cells, the pair can only be the pair it
is: observed and expected are equal by construction, whatever the biology. Such events would pull every ratio
towards 1, so they are left out of the comparison (they still count in the plain contacts). The dashboard says how
many pairs that leaves.

![Partner preference and share of contacts](../images/dash_interactions.png)

??? note "Technical: the expected contacts, exactly"

    No random shuffling is done: the expectation is computed exactly. In an event of *n* cells, of which $c_A$ are of
    type A, a given touching pair consists of an A and a B cell (A ≠ B) with probability

    $$p_{AB} = \frac{2\,c_A\,c_B}{n\,(n-1)}$$

    and of two A cells with probability $c_A (c_A - 1) / (n (n-1))$. The expected contacts of A with B are the sum over
    the events of three or more cells of (number of pairs) × $p_{AB}$, counted from both sides as the observed ones
    (an A–A pair twice). So the page and `Interactions.tsv` give identical numbers, and type pairs that never touch
    but would by chance also get a row (observed 0: avoidance).

    The sum of the expected contacts equals the sum of the observed ones: chance only redistributes the contacts over
    the types.

### What the selection in IDEAS means

The export contains the objects that were selected in IDEAS: often only multiplets, without single cells, and
sometimes only events with both cell types. The comparison with chance accounts for which cells are in each event
(it shuffles only within the event), but not for *why* an event was exported. So it answers:

> Within the exported multiplets of three or more cells, which cell types touch each other more, or less, than random
> pairing of the cells in their events would give?

It does **not** answer how often immune cells bind tumor cells in the sample: that needs the cells that never touched
(the singlets) and the events without contacts, which are not in the export. A conjugate frequency ("the percentage of
T cells in conjugates") can only be computed from an export that keeps all events, singlets included.

??? note "Technical: further limits of the null model"

    - The model treats every pair of cells in an event as equally likely to touch. In reality, large cells have more
      room for neighbours than small ones, so in an event with one large and several small cells the large cell takes
      part in more contacts than chance gives it. Comparing the same cell types between conditions is safer than
      reading one ratio on its own.
    - Small numbers make large ratios: a ratio of 2 from 2 observed and 1 expected contact means little. The tooltip
      in the dashboard and the table give the counts; look at them.
    - Artefacts take part like cells; with a stricter artefact rule they move to the cell types, and back.
    - With no preference at all, *tumor–tumor* came out at about 1.5 in tests (larger cells, more neighbours), and real
      preferences come out weaker than they are: in a test, a subset placed next to tumor cells far more often after a
      treatment rose only from 1.03 to 1.11, consistently in both donors. Compare datasets (donors, conditions) with
      each other rather than a ratio with 1.

## Comparing conditions

The **experiment** label does two things: datasets with the same label are classified together (one threshold per
marker, see [pools](classification.md#several-datasets-pools)), and the overview sums them into one partner-preference
matrix. To compare, for example, control and treated samples stained in the same way:

- give them **one** experiment label, so that "positive for A" means the same in both (or set the thresholds by hand);
- compare the conditions **per dataset**: every dataset's own table is in `interactions.tsv` (and in its own
  `Interactions.tsv`), with the observed and expected contacts; look at the same ratio in every donor, and at the
  counts behind it.

## Multiplets

`Multiplets.tsv` counts the events with *n* cells and with *n* contacts (*n* = 0, 1, 2, …), as the *multiplets* table
of ClassifyCellTypes 2.8.
