# Cell contacts

## Which cells touch

Two cells of the same event **touch** when one comes within *Membrane / interface width* of the other (empty, the
default: about 1.5 µm, i.e. 2 px at 20×, 3 px at 40×, 5 px at 60×): the same distance that defines the
[interfaces](membrane.md). So a cell with an interface always lists its partner among
the cells it touches. Cells in different events never touch.

Every cell lists the cells it touches (`Touching (labels)`) and their types (`Neighbour types`); every event counts its
pairs (`Interactions`).

??? note "Technical: the distance"

    "Within *k* pixels" means *k* steps between 4-connected pixels: the reach of the dilation that makes the
    interfaces (as CLIJ2 in the macro). Two cells with no background pixel between them are 1 step apart; with a gap of
    *g* pixels, *g* + 1. ClassifyCellTypes 2.8 used MorphoLibJ's region adjacency graph on the shrunk labels, which
    counts only cells that touch directly; here cells that almost touch also count, as they do for the interfaces.

    The default is set in µm, not in pixels, so that "touching" means the same at every magnification: the number of
    pixels nearest $1.5\,\mu\mathrm{m}$ / pixel size, at least 2. A fixed 3 px meant 3 µm at
    20× (near misses counted as contacts) and 1 µm at 60×. A width set in pixels is used as it is.

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

## Subset preference: which immune cells touch

The partner preference above shuffles all cell types, and so compares tumor cells with immune cells, which differ in
size. For the immune subtypes a cleaner question is: **when an event offers a choice of immune cells, is the one that
touches the tumor cell (or another immune cell) more often of this subtype than chance?** For this only the subtypes
are shuffled, among the immune cells of each event; the cells, their positions, the tumor cells and the contacts stay
as they are. Only events with at least two immune cells take part: with one immune cell there is no choice.

The result is observed / expected per subtype and partner (tumor cells or a subtype), with a **95 % interval**; with
fewer than 5 contacts expected there is no interval. In the dashboard (*Interactions*, tab *among immune cells*) only
squares whose interval excludes 1 are coloured. Many squares are tested at once: one square out of twenty excludes 1
by chance alone, so trust a preference that holds in every dataset (donor, replicate), not one square. The tables
are `Subset preference.tsv` per dataset and `subset preference.tsv` in the project folder.

![The subset preference in the dashboard](../images/dash_subsets.png)

??? note "Technical: the subset preference, exactly"

    **Expected.** In an event with *n* immune cells, *c<sub>S</sub>* of subtype *S*, every immune cell has subtype *S*
    with probability *c<sub>S</sub>* / *n* when the subtypes are shuffled among them. So each immune–tumor contact of the
    event is one of subtype *S* with that probability, and each immune–immune contact is one of subtypes *S* and *T*
    with probability 2 *c<sub>S</sub>* *c<sub>T</sub>* / (*n* (*n* − 1)), or *S* and *S* with probability
    *c<sub>S</sub>* (*c<sub>S</sub>* − 1) / (*n* (*n* − 1)), counted from both sides as in the contact table. The expected
    contacts are these probabilities summed over the contacts of all events with *n* ≥ 2; observed are the contacts
    counted in the same events (`contacts`). Artefacts take no part; tumor cells are never shuffled.

    **The interval.** The events are resampled with replacement (as many as there are), observed and expected are
    summed again, and observed / expected is recomputed; the 2.5th and 97.5th percentiles over 1000 resamples (Python,
    seeded, so the tables are reproducible) or 400 (the page; slightly different) are the interval. It reflects the
    sampling of events within one dataset, not the variation between donors: compare datasets for that. Below 5
    expected contacts no interval is given: a resample of a handful of contacts often contains the same contacts
    every time, which gives an interval of zero width that excludes 1 on one or two contacts.

    **Why within events.** Shuffling the subtypes over the whole dataset asks whether a subtype is over-represented
    among the immune cells that touch tumor cells. In an export of conjugates nearly every immune cell touches a tumor
    cell (95 % in the example data), so that ratio can hardly exceed about 1.05 whatever the biology. Within an event
    with a choice it can. In the example data (Subtype A+ cells were placed next to the tumor cell more often) the
    dataset-wide version gave 1.00 (0.93–1.08), the within-event version 1.09 (1.03–1.17). Effects are small in both
    because most events offer little choice: one tumor cell with two immune cells usually lets both touch it.

    **What it does not remove.** A subtype that is larger, or more often segmented as two cells, still touches more;
    subtypes are assumed to be of similar size. The selection in IDEAS applies as above.

## Comparing conditions

The **experiment** label does two things: datasets with the same label are classified together (one threshold per
marker, see [pools](classification.md#several-datasets-pools)), and the overview sums them into one partner-preference
matrix. To compare, for example, control and treated samples stained in the same way:

- give them **one** experiment label, so that "positive for A" means the same in both (or set the thresholds by hand);
- give every dataset its **condition** and **replicate** (donor): the overview then tests the change of the subset
  preference across the donors, with an interval and a correction for the number of pairs tested. See
  [Comparing conditions](comparing-conditions.md). The overview's matrix per experiment pools the conditions and
  cannot show a change;
- expect **small numbers**: few events give the immune cells a choice, so real changes of preference come out weak.
  In a test where a subset was placed next to the tumor cell more than twice as often after treatment, the ratio
  rose by 0.05–0.25, in the data's own geometry by about as much: the tool cannot recover more than the events
  contain. How many donors and events such a change needs is in
  [Comparing conditions](comparing-conditions.md#technical-calibration-by-simulation).

## Multiplets

`Multiplets.tsv` counts the events with *n* cells and with *n* contacts (*n* = 0, 1, 2, …), as the *multiplets* table
of ClassifyCellTypes 2.8.
