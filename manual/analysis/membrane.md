# Membrane enrichment

This analysis asks whether the marker of interest is enriched in the membrane of one cell type at the places where it
touches another cell type, compared with the rest of its membrane. Tick *Membrane profiles and enrichment* to switch
it on.

## The membrane profile

For every cell, the outline is traced, smoothed and **straightened**. The marker of interest is sampled in a band
across the outline, and at every position along the outline the brightest sample across the band is kept. The width
of the band is set by *Straightening width*; by default it is about 4.5 µm, which is 5 px at 20×, 9 px at 40× and
13 px at 60×. The result is a background-subtracted profile along the membrane, with one value per pixel of perimeter.

Along the same line, every position also records whether another cell lies within the interface width, and if so,
which cell. This record is the cell's **neighbor track**.

??? note "Technical: one membrane per cell"

    The Fiji macro straightened one ring for each blob of touching cells of the same type. Here every cell has its own
    membrane, so two touching tumor cells give two membranes. The outline follows the pixel edges (marching squares)
    and is smoothed with a Gaussian (σ = 2 samples), where Fiji fits a spline; it is then resampled at 1 px. At every
    pixel of the cell that lies within the interface width of another cell, the neighbor track stores the label of the
    nearest other cell. The length of every profile is stored, so a value of 0 is data and does not mark the end of the
    profile (in the macro, zeros ended it).

??? note "Technical: the straightening width, and a false enrichment it prevents"

    Where a cell is free, the outline that Cellpose draws lies a little **outside** the stained membrane (by about
    1 µm in the tests, after *Shrink cells*). Where two cells touch, it lies on the contact line. The band across the
    outline must reach the membrane in both places. If the band is too narrow, the brightest sample at the free
    membrane misses part of the membrane, while at a contact it does not. The interface then looks enriched for any
    membrane marker, without any biological cause.

    In a test at 40× (0.5 µm per pixel), a ligand on the tumor membrane was enriched only where Beta+ immune cells
    touched. The tumor cells' interface with Beta− cells, which has no enrichment, measured:

    | Straightening width | Beta− contacts (truth: none) | Beta+ contacts (truth: enriched) |
    |---|---|---|
    | 5 px (2.5 µm; the macro's value) | 1.4–1.6 SD, ratio 1.31 | 7.2–7.5 SD, ratio 2.5–2.65 |
    | 9 px (4.5 µm) | 0.35 SD, ratio 1.05 | 6.5–6.7 SD, ratio 1.95–2.0 |
    | 13 px (6.5 µm) | 0.2 SD, ratio 1.03 | 6.6–7.1 SD, ratio 1.9 |

    The default is therefore about 4.5 µm: the odd number of pixels nearest to $4.5\,\mu\mathrm{m}$ / pixel size,
    with a minimum of 5. A wider band somewhat lowers the ratio of a real enrichment, because the free membrane is now
    measured at its true brightness, but the real enrichment stays clearly separate from the false one. If you set
    the width in pixels (`membrane.straightening_width`), that value is used as it is. The comparison with the Fiji
    macro uses the macro's 5 px.

    A wider band also reaches further into a touching partner. This matters only when the partner carries the marker
    (see [Where the marker is expressed](#where-the-marker-is-expressed)).

## A pair of cell types

**Membrane of** and **Interface with** choose the pair of cell types. Each can be:

- a cell type (*Tumor*, *Immune*);
- an immune subtype as it is named in the tables (*Immune Subtype A+*, or *Immune negative* for the immune cells that
  are negative for every subtype marker);
- a marker (*Subtype A+*: every immune cell that is positive for that marker, whatever its other markers).

If a field is left empty, the cells that express the marker of interest (set in *Marker of interest expressed in*) are
measured against the other cell type, as in the macro. In the [dashboard](../guide/dashboard.md#membrane-enrichment)
you can change the pair directly.

For the measured cells, the **interface** is the part of the membrane where the neighbor track shows a partner. The
**free membrane** is the part that touches no cell at all.

## The enrichment

For every measured cell that has an interface, the analysis reports:

- **interface Mean** and **non-interface Mean**: the mean of the profile at the interface and on the free membrane,
  and the ratio of the two;
- **stdDevs above mean**: how far the interface lies above the free membrane, expressed in standard deviations of the
  free membrane:

$$\text{SD above mean} = \frac{\overline{\text{interface}} - \overline{\text{free}}}{\text{SD}(\text{free})}$$

The profiles are normalized in the same way, and aligned so that the center of the longest interface of each cell lies
at position 20. Their average, over the cells with an interface only, shows the shape of the enrichment.

![Membrane enrichment in the dashboard](../images/dash_membrane.png)

??? note "Technical: why the free membrane is the baseline"

    The macro compared the interface with the rest of the membrane, including the parts that touch other cells. A
    contact with another cell that carries the marker (for example between two tumor cells, when the marker is on tumor
    cells) places two membranes side by side and is therefore brighter. It would raise the baseline and lower the
    enrichment. Here the baseline is the membrane that touches no cell. On validation data, this raised the median
    enrichment by about 0.2 SD, while the control stayed the same.

## The control

Cells without an interface get **control windows**: *Control windows per cell* (default 10) windows at random positions
on their free membrane. Each window is as long as the median interface and is scored like an interface. The histogram
of the controls shows how much a stretch of free membrane varies by chance in this dataset.

Both histograms (in the dashboard and in the files `Histogram_interface_stdDevs_above_mean.tsv` and
`Histogram_control_stdDevs_above_mean.tsv`) start at −5 SD and share one scale. If *Histogram maximum* and
*Histogram bin width* are left empty, the scale is automatic. The maximum is the 99th percentile of the interface
values or of the control windows, whichever is higher, but at least 5 SD, and it is rounded up to a whole number of
bins. The bin width is the smallest of 0.1, 0.2, 0.25, 0.5, 1, 2.5, 5, 10, and so on, that gives at most 60 bins.
Values above the maximum are not counted. You can change the scale in the dashboard without running the analysis
again.

**The control windows are not a control for the interface itself.** Because they lie on free membrane, they cannot
show what the contact does to the measurement, independently of the marker. A contact can change the measurement in
several ways: two stained membranes lie side by side (see below), the outline runs a little outside the membrane where
the cell is free but not at a contact, and the partner's signal blurs into the interface. In a test on synthetic data
without any enrichment, the outline effect alone gave a ratio of about 1.3 and 1–1.5 SD with the macro's
straightening width, and very small p-values against the control windows. The default width now removes most of this
effect (see the technical note on the straightening width), but not the other effects. Therefore:

- **Compare pairs**, not an interface with the control windows. Compare the same measured cells with partners that do
  carry a receptor or ligand and with partners that do not (for example *Immune A+* against *Immune negative*). Or
  compare cells in which you expect the marker to be enriched with cells of the same type (set *Interface with* to
  the measured type itself). The dashboard shows all pairs at once (*All pairs* in the membrane card). In the example
  below, the tumor cells' membrane is enriched only at contacts with *Subtype A+* cells.

![All pairs of the membrane card](../images/dash_pairs.png)

- A small shift (a ratio below about 1.5, or 1–2 SD) is within what geometry alone can produce; a large effect that
  is specific to one pair is not.
- Check in the dashboard that the outlines lie on the membrane (*Events*, overlay *outlines*).
- A p-value from thousands of cells against thousands of windows says nothing about biological replicates. Compare
  datasets (donors, replicates) instead. A change between conditions is tested across the donors in the overview; see
  [Comparing conditions](comparing-conditions.md).

??? note "Technical: comparing pairs, and why there is no reference channel"

    **The pair table.** For every selection under *Membrane of* (rows) and every cell type or subtype under
    *Interface with* (columns), the dashboard computes the median *SD above mean* of the measured cells that have an
    interface with that partner, and the number of such cells. The median is computed exactly as for the chosen pair,
    but without the control windows. A pair in which the measured cells or the partners do not carry the marker is the
    geometry-matched control for a pair in which they do: the same cells, the same kind of contact and the same
    segmentation. Contacts between cells of the same type are the control for two stained membranes side by side:
    their enrichment shows what the doubling alone produces. Each entry carries the marker check of its pair: "!" when
    the partners carry part of the marker or as much as the measured cells, "!!" when they carry more, and "–" when the
    measured cells hardly carry it.

    **No reference channel.** A second channel, measured like the marker of interest but not expected to be enriched,
    seemed a natural control. The cells' own tracker dye was tested as such a reference in the 40× test described
    above. Because the dye is cytoplasmic, it fills the band at a contact but not at the free edge, and it showed 4.8 SD
    (at 5 px) and 2.7 SD (at 9 px) of "enrichment". A reference must be stained in the same way as the marker (a
    membrane stain on the same cells), which a panel rarely includes. The pair comparison needs no extra channel.

??? note "Technical: the control windows"

    The windows are drawn on the free membrane only, which is read as one closed line, and each window is scored
    against the rest of that free membrane. The macro used one window of 25 px per membrane; many windows of the
    typical interface length give a better estimate of the null distribution. The windows are random and not seeded,
    so the control differs slightly from run to run (the dashboard uses a fixed seed). All windows are listed in
    `Control windows stdDevs above mean.tsv`.

## Where the marker is expressed

**The marker of interest must be on the measured cells, not on their partners.** At an interface, the membranes of the
two cells lie against each other. If the partner also carries the marker, its membrane is measured together with that
of the measured cell, and the interface looks enriched even when both cells are stained evenly.

The analysis checks this. It compares the median free-membrane intensity of the partners with that of the measured
cells. The result is reported in the log, in `run.json`, in the overview (*MOI partners / measured*) and as a badge in
the dashboard:

| Partners / measured | Badge | Meaning |
|---|---|---|
| below 25 % | ✓ green | The marker is mainly on the measured cells. |
| 25–80 % | ! amber | The partners carry part of the marker, so part of the interface signal may come from them. |
| 80–125 % | ! amber | The marker is on both cell types at a similar level. Two stained membranes side by side look like enrichment; compare with contacts between cells of the same type. |
| 125 % or more | ! red | The marker is brighter on the partners. Check *Marker of interest expressed in* and the pair. |
| (any) | – gray | The measured cells hardly carry the marker (below 25 % of the brightest cell type), so their enrichment is not informative. |

??? note "Technical: the marker check, exactly"

    Every cell gets a **free-membrane intensity**: the mean of its straightened profile where its neighbor track shows
    no cell (background-subtracted). Without the membrane analysis, the mean of the whole cell is used instead.
    *Measured* and *partners* are the medians of this intensity over the cells of each selection. The **brightest cell
    type** is the highest median over the cell types and immune subtypes, excluding artefacts. The rules are applied
    in this order:

    1. partners / measured ≥ 125 %: brighter on the partners (red);
    2. measured below 25 % of the brightest cell type: not informative (gray);
    3. 80–125 %: on both (amber);
    4. 25–80 %: the partners carry part (amber);
    5. below 25 %: mainly on the measured cells (green).

    These limits are rules of thumb. At 25 %, the partner adds a quarter of the measured cell's membrane signal to
    every interface pixel it reaches. At 80–125 %, an interface between two equally stained membranes looks like a
    doubling, whatever the biology.

If the badge is amber or red, you can choose the pair the other way round (if the marker is on the partners), compare
the enrichment with that of a pair whose partners do not carry the marker, or treat the result as an upper limit.

**Look at the subsets, on both sides.** The default pair (all cells of the type that expresses the marker, with all
cells of the other type) can hide an effect that only one subset shows. In one test, a receptor that was strongly
polarized on *A+* immune cells gave a median of only 0.2 SD over all immune cells, and a ligand that was enriched only
where *B+* cells touched gave 0.7 SD with all immune cells. The *All pairs* table in the dashboard and
`Membrane pairs.tsv` give every combination of *Membrane of* and *Interface with* at once. Read them before you conclude
that there is no enrichment.

??? note "Technical: inside only, and why it is no cure"

    *Straighten inside the measured cell only* takes the brightest sample only from the part of the band that lies
    inside the measured cell, so that the band does not reach into the partner. This is a **sensitivity test, not a
    separation**. Touching membranes are nanometers apart, far below the resolution of the ImageStream, so the blur
    mixes both membranes in every pixel near the interface, whichever pixels are sampled. If the enrichment drops
    clearly with *inside only*, part of it comes from the partner side. If it stays the same, that is reassuring but
    not proof.

## The per-event table

Independently of the profiles, every event gets the mean of the marker of interest at the interface and on the rest of
the membrane, for each cell type. These values are in `All cell combinations.tsv`, and for the events with one tumor
cell and one immune cell also in `Single Tumor Single Immune cell.tsv`, as in the macro. They are measured in the
membrane band of each cell type (*Membrane / interface width*), not in the straightened profiles.

This is the macro's measure, and it has none of the safeguards described above. It takes all partners of the other
type together, has no control, and does not follow the chosen membrane pair. A partner that carries the marker, or two
stained membranes side by side, raise it as much as a real enrichment does (in a test without any enrichment it was
1.27). The overview column *median … membrane ratio (doublets)* and the dashboard card *Marker of interest at the
interface* show this measure. Read them together with the marker check, and base your conclusions on the pair
comparisons.
