# Membrane enrichment

The question: is the marker of interest enriched in the membrane of one cell type where it touches another, compared
with the rest of its membrane? *Membrane profiles and enrichment* switches this analysis on.

## The membrane profile

For every cell, the outline is traced, smoothed, and **straightened**: the marker of interest is sampled in a band
across the outline (*Straightening width*, 5 px), and the brightest sample across the band is kept at every position.
This gives a profile along the membrane, one value per pixel of perimeter, background-subtracted.

Along the same line, every position records whether another cell is within the interface width, and which one: the
cell's **neighbour track**.

??? note "Technical: one membrane per cell"

    The Fiji macro straightened one ring per blob of touching cells of the same type; here every cell has its own
    membrane, so two touching tumor cells are two membranes. The outline runs along the pixel edges (marching squares),
    smoothed with a Gaussian (σ = 2 samples) where Fiji fits a spline, and resampled at 1 px. The neighbour track
    takes, at every own pixel within the interface width of another cell, the label of the nearest other cell. The
    length of every profile is kept; a value of 0 is data, not the end of the profile (in the macro, zeros ended it).

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
In tests on synthetic data without any enrichment, such effects alone gave a ratio of about 1.3 and 1–1.5 SD, with
very small p-values against the control windows. So:

- **Compare pairs**, not an interface with the control windows: the same measured cells with partners that do carry a
  receptor or ligand and partners that do not (for example *Immune A+* against *Immune negative*), or cells whose
  marker you expect to be enriched against cells of the same type (*Interface with* the measured type itself).
- A small shift (ratio below about 1.5, or 1–2 SD) is within what geometry alone can give; a large, pair-specific
  effect is not.
- Check in the dashboard that the outlines lie on the membrane (*Events*, overlay *outlines*).
- The p-value of thousands of cells against thousands of windows says nothing about biological replicates: compare
  datasets (donors, replicates) instead.

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

If the badge is amber or red, choose the pair the other way round (if the marker is on the partners), compare the
enrichment with that of a pair whose partners do not carry the marker, or treat the result as an upper bound.

**Look at the subsets.** The default pair (all cells of the type that expresses the marker) can hide an effect that
only one subset has: in a test, a receptor polarized strongly on *A+* immune cells gave a median of 0.2 SD over all
immune cells. Choose each subset under *Membrane of* in the dashboard before concluding that there is no
enrichment.

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
