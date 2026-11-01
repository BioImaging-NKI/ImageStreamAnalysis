# Membrane enrichment

The question: is the marker of interest enriched in the membrane of one cell type where it touches another, compared
with the rest of its membrane? *Membrane profiles and enrichment* switches this analysis on.

## The membrane profile

For every cell, the outline is traced, smoothed, and **straightened**: the marker of interest is sampled in a band
across the outline (*Straightening width*: about 4.5 µm, i.e. 5 px at 20×, 9 px at 40×, 13 px at 60×), and the
brightest sample across the band is kept at every position.
This gives a profile along the membrane, one value per pixel of perimeter, background-subtracted.

Along the same line, every position records whether another cell is within the interface width, and which one: the
cell's **neighbour track**.

??? note "Technical: one membrane per cell"

    The Fiji macro straightened one ring per blob of touching cells of the same type; here every cell has its own
    membrane, so two touching tumor cells are two membranes. The outline runs along the pixel edges (marching squares),
    smoothed with a Gaussian (σ = 2 samples) where Fiji fits a spline, and resampled at 1 px. The neighbour track
    takes, at every own pixel within the interface width of another cell, the label of the nearest other cell. The
    length of every profile is kept; a value of 0 is data, not the end of the profile (in the macro, zeros ended it).

??? note "Technical: the straightening width, and a false enrichment it prevents"

    The outline that Cellpose draws lies a little **outside** the stained membrane where a cell is free (by about
    1 µm in the tests, after *Shrink cells*), but on the contact line where two cells touch. The band across the
    outline must reach the membrane in both places. With a band that is too narrow, the brightest sample at the free
    membrane misses part of the membrane, while at a contact it does not: the interface then looks enriched for any
    membrane marker, with no biology.

    In a test at 40× (0.5 µm per pixel; a ligand on the tumor membrane, enriched only where Beta+ immune cells touch),
    the tumor cells' interface with Beta− cells, which has no enrichment, measured:

    | Straightening width | Beta− contacts (truth: none) | Beta+ contacts (truth: enriched) |
    |---|---|---|
    | 5 px (2.5 µm; the macro's value) | 1.4–1.6 SD, ratio 1.31 | 7.2–7.5 SD, ratio 2.5–2.65 |
    | 9 px (4.5 µm) | 0.35 SD, ratio 1.05 | 6.5–6.7 SD, ratio 1.95–2.0 |
    | 13 px (6.5 µm) | 0.2 SD, ratio 1.03 | 6.6–7.1 SD, ratio 1.9 |

    So the default is about 4.5 µm: the odd number of pixels nearest $4.5\,\mu\mathrm{m}$ / pixel size, at least 5.
    A wider band lowers the ratio of a real enrichment somewhat (the free membrane is now measured at its true
    brightness) but keeps it clearly apart. A width set in pixels (`membrane.straightening_width`) is used as it is;
    the comparison with the Fiji macro uses its 5 px.

    A wider band also reaches further into a touching partner. That matters only when the partner carries the marker
    (see [Where the marker is expressed](#where-the-marker-is-expressed)).

## A pair of cell types

**Membrane of** and **Interface with** choose the pair. Each can be a cell type (*Tumor*, *Immune*), an immune subtype
as named in the tables (*Immune Subtype A+*, or *Immune negative* for the immune cells negative for every subtype
marker), or a marker (*Subtype A+*: every immune cell positive for it, whatever
its other markers). Empty: the cells that express the marker of interest (*Marker of interest expressed in*) against
the other type, as in the macro. In the [dashboard](../guide/dashboard.md#membrane-enrichment) the pair can be changed
on the spot.

For the measured cells, the **interface** is the part of the membrane whose neighbour track shows a partner; the
**free membrane** is the part that touches no cell at all.

## The enrichment

Per measured cell with an interface:

- **interface Mean** and **non-interface Mean**: the mean of the profile at the interface and on the free membrane,
  and their ratio;
- **stdDevs above mean**: how far the interface lies above the free membrane, in standard deviations of the free
  membrane:

$$\text{SD above mean} = \frac{\overline{\text{interface}} - \overline{\text{free}}}{\text{SD}(\text{free})}$$

The profiles are normalized the same way and aligned so that the centre of the longest interface lies at position 20;
their average (cells with an interface only) shows the shape of the enrichment.

![Membrane enrichment in the dashboard](../images/dash_membrane.png)

??? note "Technical: why the free membrane is the baseline"

    The macro compared the interface with the rest of the membrane, including the parts that touch other cells. A
    contact with another cell that carries the marker (for example two tumor cells, when the marker is on tumor cells)
    places two membranes side by side and is brighter; it would raise the baseline and lower the enrichment. Here the
    baseline is the membrane that touches no cell. On validation data this raised the median enrichment by about
    0.2 SD; the control stayed the same.

## The control

Cells without an interface get **control windows**: *Control windows per cell* (10) windows at random positions on
their free membrane, as long as the median interface, scored like an interface. The histogram of the controls shows
how much a stretch of free membrane varies by chance in this dataset.

**The control windows are not a control for the interface itself.** They lie on free membrane, so they cannot show
what the contact does to the measurement regardless of the marker: two stained membranes side by side (see below), an
outline that runs a little outside the membrane where the cell is free but not at a contact, or blur from the partner.
In a test on synthetic data without any enrichment, the outline effect alone gave a ratio of about 1.3 and
1–1.5 SD with the macro's straightening width, and very small p-values against the control windows; the default
width now removes most of it (see the technical note on the straightening width), but not the other effects. So:

- **Compare pairs**, not an interface with the control windows: the same measured cells with partners that do carry a
  receptor or ligand and partners that do not (for example *Immune A+* against *Immune negative*), or cells whose
  marker you expect to be enriched against cells of the same type (*Interface with* the measured type itself). The
  dashboard shows all pairs at once (*All pairs* in the membrane card): in the example below, the tumor cells'
  membrane is enriched only at contacts with *Subtype A+* cells.

![All pairs of the membrane card](../images/dash_pairs.png)

- A small shift (ratio below about 1.5, or 1–2 SD) is within what geometry alone can give; a large, pair-specific
  effect is not.
- Check in the dashboard that the outlines lie on the membrane (*Events*, overlay *outlines*).
- The p-value of thousands of cells against thousands of windows says nothing about biological replicates: compare
  datasets (donors, replicates) instead.

??? note "Technical: comparing pairs, and why there is no reference channel"

    **The pair table.** For every selection under *Membrane of* (rows) and every cell type or subtype under
    *Interface with* (columns), the dashboard computes the median *SD above mean* of the measured cells that have an
    interface with that partner, exactly as for the chosen pair but without the control windows, and the number of
    such cells. A pair whose measured cells or partners do not carry the marker is the geometry-matched control of a
    pair that does: same cells, same kind of contact, same segmentation. Contacts between cells of the same type are
    the control for two stained membranes side by side: their enrichment is what doubling alone gives. Each entry
    carries the marker check of its pair (! partners carry part of the marker or as much, !! more, – the measured
    cells hardly carry it).

    **No reference channel.** A second channel measured like the marker of interest, one that should not be enriched,
    looked like a natural control. The cells' own tracker dye was tested as such a reference in the 40× test above:
    being cytoplasmic, it fills the band at a contact but not at the free edge, and showed 4.8 SD (5 px) and 2.7 SD
    (9 px) of "enrichment". A reference must be stained like the marker (a membrane stain on the same cells), which a
    panel rarely has; the pair comparison needs no extra channel.

??? note "Technical: the control windows"

    The windows are drawn on the free membrane only, read as one closed line, and scored against the rest of that
    free membrane. The macro used one window of 25 px per membrane; many windows of the typical interface length
    estimate the null distribution better. The windows are random and not seeded, so the control differs slightly
    from run to run (the dashboard uses a fixed seed). All windows are in `Control windows stdDevs above mean.tsv`.

## Where the marker is expressed

**The marker of interest must be on the measured cells, not on their partners.** At an interface, the membranes of the
two cells lie against each other. If the partner carries the marker too, its membrane is measured along with the
measured cell's, and the interface looks enriched even when both cells are stained evenly.

The analysis checks this. It compares the median free-membrane intensity of the partners with that of the measured
cells, and reports it in the log, `run.json`, the overview (*MOI partners / measured*) and as a badge in the dashboard:

| Partners / measured | Badge | Meaning |
|---|---|---|
| below 25 % | ✓ green | the marker is mainly on the measured cells |
| 25–80 % | ! amber | the partners carry part of it: part of the interface signal may be theirs |
| 80–125 % | ! amber | the marker is on both at a similar level: two stained membranes side by side look like enrichment; compare with contacts between cells of the same type |
| 125 % or more | ! red | the marker is brighter on the partners: check *Marker of interest expressed in* and the pair |
| (any) | – grey | the measured cells hardly carry the marker (below 25 % of the brightest cell type): their enrichment is not informative |

??? note "Technical: the marker check, exactly"

    For every cell its **free-membrane intensity**: the mean of its straightened profile where its neighbour track
    shows no cell (background-subtracted; without the membrane analysis the mean of the whole cell). *Measured* and
    *partners* are the medians over the cells of each selection; the **brightest cell type** is the highest median
    over the cell types and immune subtypes (artefacts excluded). In this order: partners / measured ≥ 125 %:
    brighter on the partners (red); measured below 25 % of the brightest type: not informative (grey); 80–125 %: on
    both (amber); 25–80 %: the partners carry part (amber); below 25 %: mainly on the measured cells (green). The
    bands are rules of thumb: at 25 % the partner adds a quarter of the measured cell's membrane signal to every
    interface pixel it reaches; at 80–125 % an interface between two equally stained membranes looks like a doubling
    whatever the biology.

If the badge is amber or red, choose the pair the other way round (if the marker is on the partners), compare the
enrichment with that of a pair whose partners do not carry the marker, or treat the result as an upper bound.

**Look at the subsets, on both sides.** The default pair (all cells of the type that expresses the marker, with all
cells of the other type) can hide an effect that only one subset has: in a test, a receptor polarized strongly on *A+*
immune cells gave a median of 0.2 SD over all immune cells, and a ligand enriched only where *B+* cells touch gave
0.7 SD with all immune cells. The *All pairs* table of the dashboard and `Membrane pairs.tsv` give every combination of
*Membrane of* and *Interface with* at once; read them before concluding that there is no enrichment.

??? note "Technical: inside only, and why it is no cure"

    *Straighten inside the measured cell only* takes the brightest sample only from the part of the band that lies
    inside the measured cell, so the band does not reach into the partner. This is a **sensitivity test, not a
    separation**: touching membranes are nanometres apart, far below the resolution of the ImageStream, so the blur
    mixes both membranes in every pixel near the interface, whichever pixels are sampled. If the enrichment drops
    clearly with *inside only*, part of it comes from the partner side; if it stays, that is reassuring but not proof.

## The per-event table

Independently of the profiles, every event gets the mean of the marker of interest at the interface and on the rest
of the membrane of each cell type (`All cell combinations.tsv`; for the events with one tumor and one immune cell,
`Single Tumor Single Immune cell.tsv`), as in the macro. These use the membrane band of each cell type
(*Membrane / interface width*) rather than the straightened profiles.

This is the macro's measure, and it has none of the safeguards above: it takes all partners of the other type
together, has no control, and does not follow the membrane pair. A partner that carries the marker, or two stained
membranes side by side, raise it as much as a real enrichment does (in a test without any enrichment it was 1.27).
The overview's *median … membrane ratio (doublets)* and the dashboard card *Marker of interest at the interface*
show it; read them together with the marker check, and use the pair comparisons for conclusions.
