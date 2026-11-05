# Membrane enrichment

This analysis asks a simple question: **is there more of a marker in a cell's membrane where it touches another
cell than in the rest of its membrane?** For example: does a ligand on tumor cells gather where a T cell touches?
Tick *Membrane profiles and enrichment* to switch it on. With several markers of interest, every marker is analyzed
in the same way (see [Several markers of interest](#several-markers-of-interest)).

The answer is less simple than the question. Images of touching cells contain several effects that look like
enrichment but are not: two stained membranes lying side by side, a partner that carries the marker too, and an
outline that does not lie exactly on the membrane. This page explains each step, what can go wrong, and what the
program does about it. The main text is meant to be read from top to bottom; the folded *Technical* notes give the
details for those who want them.

**In short**, for every cell the program:

1. traces the cell's outline and finds where the stained membrane lies relative to that outline;
2. reads the marker along the membrane, as a line of values from one end of the membrane to the other (the
   **membrane profile**);
3. notes, for every point of that line, whether another cell touches there, and which one;
4. compares the part that touches a partner (the **interface**) with the part that touches no cell at all (the
   **free membrane**).

The result is a score per cell, and a median per pair of cell types. Random stretches of free membrane serve as a
control, and contacts between cells of the same kind serve as a reference for what touching alone does.

## The membrane profile

Imagine laying a ribbon along the edge of the cell, and reading the marker under the ribbon at every step along the
edge. That is what the program does: it traces the cell's outline, smooths it, and takes a **band** of pixels across
the membrane at every point along it. Laid out straight, the band becomes a strip; this is called *straightening*.
At every point along the membrane, the brightest pixel across the band is kept. The result is the membrane profile:
one value per pixel of perimeter, with the background subtracted.

The width of the band is *Straightening width*. By default it is about 4.5 µm (5 px at 20×, 9 px at 40×, 13 px at
60×). It has to be wide enough to contain the membrane everywhere, also where the outline does not lie exactly on it
(see the next section).

Along the same line, the program also records, for every point, whether another cell lies within *Membrane /
interface width* (by default about 1.5 µm), and if so which cell. This record is the cell's **neighbor track**. It
decides which points are interface and which are free membrane.

??? note "Technical: one membrane per cell, and how the outline is traced"

    The Fiji macro straightened one ring for each blob of touching cells of the same type. Here every cell has its own
    membrane, so two touching tumor cells give two membranes. The outline follows the pixel edges of the cell's mask
    (marching squares) and is smoothed with a Gaussian (σ = 2 samples), where Fiji fits a spline; it is then resampled
    at 1 px, so a profile has one value per pixel of perimeter. For the neighbor track, every pixel of the cell that
    lies within the interface width of another cell gets the label of the nearest other cell; the track is read along
    the outline. The length of every profile is stored, so a value of 0 is data and does not mark the end of the
    profile (in the macro, zeros ended it).

## Where the band is placed

### The problem: the outline is not the membrane

The outline comes from the cell's mask, which Cellpose draws. It is close to the membrane, but not on it, and **how
far off it is depends on where you look**:

- **Where a cell is free**, the mask usually reaches a little beyond the stained membrane, because the blurred edge
  of the cell is included. The membrane then lies a little *inside* the outline, typically less than a micrometre
  (about 1 µm in the synthetic tests).
- **Where two cells touch**, the outline is the border between the two masks. Cellpose has to put that border
  somewhere in the blurred contact, and it need not be exactly on the membrane either. On real data, it sometimes lay
  a little into one of the two cells.

So the distance between outline and membrane differs between the free membrane and the interface, and between cell
types. If the band is centred on the outline, it can reach more of the membrane at the interface than at the free
membrane. The interface then looks brighter than the free membrane even if the marker is spread evenly. This is a
**false enrichment**: it comes from the geometry, not from the biology.

### What the program does: find the membrane first

Before it measures anything, the program finds out where the membrane really lies. For each cell type, it averages
the marker across the outline (3 µm to each side) over all its cells: once along the free membrane, and once along
the contacts with each other cell type. The peak of each average is where the membrane lies. The band is then
centred on that peak: at the free membrane on the free peak, at a contact on the peak for that kind of contact.
Where a contact begins or ends, the band moves over smoothly.

You can see the result in the dashboard: the membrane card has a folded section *Across the membrane*, which shows
these averages and where the band lies; and in the *Cell* card, *membrane band* (or the space bar) draws the band of
a cell over its images.

*Band centred on* has three choices:

| Choice | The band lies | Use it |
|---|---|---|
| *membrane* (default) | where the membrane lies, for each cell type (one place for all cells of a type, at the free membrane and at each kind of contact) | in general |
| *cell* | as *membrane*, but each cell's free membrane on its own peak | when the cells differ a lot in how far the outline lies from the membrane |
| *outline* | on the outline, as the Fiji macro did | to compare with the macro, or to see how much the placement matters |

The *Straightening width* still matters, but much less. Centring removes the systematic difference between free
membrane and interface; the width takes up what remains (cells differ from each other). Keep the default width
unless you have a reason to change it.

??? note "Technical: how the membrane is found"

    For every measured cell, every marker of interest is sampled along lines perpendicular to the outline, from 3 µm
    inside to 3 µm outside, in steps of a quarter pixel. These samples are averaged along the points of the cell that
    are free membrane (if there are at least 5), and separately along the points that touch each cell type. The cell
    averages are then averaged over all cells of the same type. The peak of such an average (refined with a parabola
    through the three highest points) is the distance from the outline to the membrane for that cell type and place.
    A peak at the very end of the 3 µm window is not accepted.

    **Which marker places the band.** The membrane of a cell type is located with a marker of interest that this cell
    type carries: its free membrane must be at least half as bright as that of the brightest cell type (measured on the
    same averages). A marker that a cell type does not carry shows no membrane of that cell, or signal somewhere else in
    it (on real data, one marker showed a sharp peak well inside the other cell type, probably signal inside the cell).
    With several markers carried by the same type, the marker most specific to that type is used, and if two are
    equally specific, the one with the clearest peak. Its place then holds for all markers, because the membrane is
    where it is, whatever the marker.

    **When the band stays on the outline.** If no marker of interest is carried by a cell type, if fewer than 20 cells
    contribute to an average, or if the average has no clear peak (its peak must stand out by at least 30 % of its
    range above the ends of the window), the band stays on the outline there. The log says so for every cell type and
    place.

    The neighbor track is always read along the outline itself, so moving the band does not change which points are
    interface and which are free membrane.

??? note "Technical: *cell* — a place per cell, and why not at the contacts"

    With *cell*, every cell's free membrane is placed on the peak of its own average (from the same marker as its
    type), and its contacts are moved by the same amount, plus its type's difference between contact and free
    membrane. Example: if tumor cells have their free membrane on average 0.9 µm inside the outline and their contacts
    with immune cells 0.3 µm outside it, a tumor cell whose own free membrane lies 1.2 µm inside gets its contacts
    with immune cells at 0.0 µm.

    Contacts are never placed per cell, for two reasons. They are short, so the peak of a single contact is noisy. And
    letting the marker choose where to sample exactly the place that is tested would push the interface up: the
    brightest spot nearby always looks a little enriched. Cells with fewer than 20 points of free membrane, or without
    a clear peak of their own, use the type's place. The log says how many cells had a place of their own.

    Is the per-cell place real, or noise? On real datasets, the place found on one half of a cell's free membrane
    agreed well with the place found on the other half, and cells differed from each other by about as much as the
    type's average distance. *cell* follows that variation.

??? note "Technical: what was seen on real data"

    On real data, the free membrane of both cell types lay on average somewhat inside the outline (less than a
    micrometre), and single cells varied around that by about as much. At contacts, the two cell types differed. For one
    type, the membrane lay just as far inside the outline at contacts as at the free membrane, so a band on the
    outline treated both alike. For the other type, the membrane at contacts lay on or slightly beyond the outline,
    while the free membrane lay inside it. A band on the outline then reached the two places differently, and how
    that shifted the score depended on the band width: with a narrow band it could raise the score of some pairs and
    lower that of others. Centring the band on the membrane made the scores depend much less on the width.

    Contacts between cells of the same type, where no enrichment is expected, came out lower once the band was
    centred, which is what removing a false enrichment should do.

??? note "Technical: the straightening width with the band on the outline"

    In a synthetic test at 40× (0.5 µm per pixel), with the band on the outline, a ligand on the tumor membrane was
    enriched only where Beta+ immune cells touched. The tumor cells' interface with Beta− cells, which has no
    enrichment, measured:

    | Straightening width | Beta− contacts (truth: none) | Beta+ contacts (truth: enriched) |
    |---|---|---|
    | 5 px (2.5 µm; the macro's value) | 1.4–1.6 SD, ratio 1.31 | 7.2–7.5 SD, ratio 2.5–2.65 |
    | 9 px (4.5 µm) | 0.35 SD, ratio 1.05 | 6.5–6.7 SD, ratio 1.95–2.0 |
    | 13 px (6.5 µm) | 0.2 SD, ratio 1.03 | 6.6–7.1 SD, ratio 1.9 |

    This is why the default width is about 4.5 µm: the odd number of pixels nearest to 4.5 µm divided by the pixel
    size, with a minimum of 5. A wider band somewhat lowers the ratio of a real enrichment, because the free membrane
    is now measured at its true brightness, but the real enrichment stays clearly apart from the false one. If you set
    the width in pixels (`membrane.straightening_width`), that value is used as it is. The comparison with the Fiji
    macro uses the macro's 5 px with the band on the outline.

    With the band centred on the membrane, the width matters much less: in a synthetic test with an evenly stained
    membrane that lay 2 px inside a free outline, a 3 px band on the outline made the interface look about 30 %
    brighter than the free membrane; centred on the membrane, a 3 px and a 9 px band both showed no enrichment.

    A wider band also reaches further into a touching partner. This matters when the partner carries the marker
    (see [Where the marker is expressed](#where-the-marker-is-expressed)).

## A pair of cell types

**Membrane of** and **Interface with** choose the pair of cell types: whose membrane is measured, and which
partners make up the interface. Each can be:

- a cell type (*Tumor*, *Immune*);
- an immune subtype as it is named in the tables (*Immune Subtype A+*, or *Immune negative* for the immune cells that
  are negative for every subtype marker);
- a marker (*Subtype A+*: every immune cell that is positive for that marker, whatever its other markers).

Anything else (for example *Immune Subtype A* without the +) stops the run before it starts, with a message that
suggests what was probably meant.

If *Membrane of* is left empty, the program measures the cell type that **carries the marker**: the tumor cells or
the immune cells, whichever has the brighter membrane. This is decided from the data, for each marker of interest
separately. If *Interface with* is left empty, the partners are the other cell type. Every marker of interest has a
pair of its own, and each side is set on its own. In the [dashboard](../guide/dashboard.md#membrane-enrichment) you
can change the pair directly.

The log says which cell type carries each marker, and how clearly. In the example (Donor 1): "Marker: carried by the
tumor cells (...): membranes of Tumor with Immune". The same is recorded in `run.json` (*markers of interest*), in
the column *marker on* of `Membrane pairs.tsv`, and in the overview column *MOI on* of each marker. **If the two cell
types carry the marker at a similar level** (the weaker one at 80 % of the stronger or more), the choice is made by a
small margin, and the log warns. Then set the pair yourself, so that it does not depend on that margin.

??? note "Technical: how the carrying cell type is decided"

    With the membrane analysis, the measure is the height of the free membrane in the marker averaged across the
    outline (the same averages that place the band): the peak of that average above its ends, for the tumor cells and
    for the immune cells. The type with the higher peak carries the marker. This measure does not depend on the band's
    width or place. (Until November 2026 the band's own free-membrane values were used; when a marker was on both cell
    types at a similar level, the choice could then flip when only the band width changed.) Without the membrane
    analysis, the median whole-cell intensity is used instead.

    Before November 2026, the measured type was a setting (*Marker of interest expressed in*). A settings file or
    project that still has it (`channels.cell_of_interest`) loads without an error; the setting is ignored, and the
    log says so.

For the measured cells, the **interface** is the part of the membrane where the neighbor track shows a partner. The
**free membrane** is the part that touches no cell at all. Contacts with cells that are not partners belong to
neither: they are left out of the comparison.

## The score

For every measured cell that has an interface, the analysis compares the mean of the profile at the interface with
the mean of the free membrane (both after the background is subtracted). The main score is the **log2 ratio**:

$$\text{log2 ratio} = \log_2 \frac{\overline{\text{interface}}}{\overline{\text{free}}}$$

It reads as a fold change: **0** means the interface is as bright as the free membrane, **1** twice as bright, **2**
four times, and **−1** half as bright. A log2 ratio of 0.5 is about 40 % more marker at the contact.

The tables also keep the macro's score, **SD above mean**: how far the interface lies above the free membrane, in
standard deviations of the free membrane:

$$\text{SD above mean} = \frac{\overline{\text{interface}} - \overline{\text{free}}}{\text{SD}(\text{free})}$$

The profiles are aligned so that the center of the longest interface of each cell lies at position 20. Their
average, over the cells with an interface, shows the shape of the enrichment.

![Membrane enrichment in the dashboard](../images/dash_membrane.png)

??? note "Technical: why the log2 ratio, and not SD above mean"

    Both scores compare the interface with the cell's own free membrane. They differ in what they divide by: the
    ratio divides by the mean of the free membrane, the SD score by its spread.

    - **The spread mixes in other things.** The standard deviation of the free membrane grows with patchy staining
      (clusters, microvilli), with noise, and with where the band lies. A cell with an evenly stained membrane gets a
      large SD score for a small enrichment; a patchy cell gets a small score for a large one. If a treatment makes
      the marker more clustered without bringing more of it to the contact, the SD score changes and the ratio does
      not.
    - **Brightness.** Dim cells are relatively noisier, so their SD scores shrink, and bright cells get large ones.
      The ratio does not depend on how much marker a cell carries.
    - **Stability.** A mean from a short free membrane is estimated much better than a standard deviation, so the
      ratio has fewer extreme values.
    - **Statistics.** log2 ratios are symmetric (twice as bright = +1, half = −1) and close to a normal distribution,
      which suits the tests of the [comparison of conditions](comparing-conditions.md).
    - **Meaning.** "1.4 times as bright" can be judged biologically and compared with other studies; "2 SD" depends
      on the cell.

    **Its weak point is the background.** A ratio needs a correct zero: an offset that is subtracted too little or too
    much biases it, while the SD score does not care. And a cell whose free membrane is about as dark as the background
    has no meaningful ratio: the program leaves the log2 ratio empty when either mean is 0 or below. So the log2 ratio
    is reliable for the cell type that carries the marker, and much less so for a type that hardly carries it (the
    marker check marks such pairs with "–"). The SD score does say whether the interface stands out from the cell's
    own fluctuations; the control windows answer that for both scores.

??? note "Technical: why the free membrane is the baseline"

    The macro compared the interface with the rest of the membrane, including the parts that touch other cells. A
    contact with another cell that carries the marker (for example between two tumor cells, when the marker is on
    tumor cells) places two membranes side by side and is therefore brighter. It would raise the baseline and lower
    the enrichment. Here the baseline is the membrane that touches no cell. On validation data, this raised the median
    enrichment by about 0.2 SD, while the control stayed the same.

### Ways to score an interface

The dashboard's *Measure* menu offers five scores. The tables of the output folder contain all of them.

| Measure | What it compares | When to use it |
|---|---|---|
| **log2 ratio to own membrane** (default) | the interface with the cell's own free membrane, as a fold change | when the partners carry little or none of the marker (the marker check is green) |
| **log2 ratio above same-type contacts** | the same, minus the median log2 ratio of the measured cells at contacts with cells of their own kind | when the partners carry the marker too: it asks whether the contact with *this* partner is brighter than a contact between two cells of the measured kind |
| **SD above own membrane** | as the first, in standard deviations of the free membrane (the macro's score) | to compare with results of the Fiji macro |
| **SD above same-type contacts** | as the second, in standard deviations | as above |
| **SD above both membranes** (experimental) | the interface with the own free membrane plus the partner's free membrane | to explore; read the note below before relying on it |

Why the same-type versions? At an interface, the membranes of the two cells lie so close together that the
microscope cannot separate them. If the partner carries the marker too, the interface contains two stained
membranes, and it looks brighter even without any enrichment. A score against the own membrane cannot tell this apart
from a real enrichment.

The same-type versions use the data themselves as the reference. A contact between two tumor cells also places two
stained membranes side by side, with the same kind of geometry, but with no partner-specific biology. If a
tumor–immune contact scores higher than a tumor–tumor contact, that difference is what the immune partner adds. A
score of 0 means "as bright as a contact between two cells of the measured kind"; in log2, a value of 0.5 means 40 %
brighter than such a contact.

??? note "Technical: the scores above same-type contacts, exactly"

    For the measured cells of a pair (*Membrane of*), the program takes the score of each of them at its contacts with
    cells of the same selection: tumor cells with tumor cells, *Immune Subtype A+* cells with *Immune Subtype A+*
    cells, and so on. The median of these scores is the **reference**; it needs at least 10 cells with such a contact
    (otherwise the measure is empty). Each cell's score minus the reference is its score above same-type contacts. In
    log2 this is a ratio of ratios: how much brighter the contact with the partner is than a same-type contact. In the
    files: the columns `log2 ratio above same-type contacts` and `stdDevs above same-type contacts` of `Membrane
    intensity analysis - <marker>.tsv`, and `median log2 ratio above same-type contacts` and `median above same-type
    contacts` in `Membrane pairs.tsv`. In the dashboard, the histograms, the heatmap order, the *All pairs* table, the
    cell card and the *Cells* table follow the measure; the average profile stays that of the own membrane, with a
    dashed line at the reference.

    **When it is not a fair reference.** It assumes that a contact with the partner is, apart from the biology,
    like a contact between two measured cells. That holds best when the partner carries about as much of the marker as
    the measured cells (the marker check says *on both cells*). If the partner carries much less, the same-type
    contacts contain more marker than a fair reference would, and the measure is too strict: it can hide an
    enrichment. If the partner carries more, the reference is too lenient. Cells of different sizes also touch in
    different ways. Read it together with the marker check, and with *SD above own membrane*.

??? note "Technical: SD above both membranes (experimental), and why it fails its own check"

    The idea: if the intensities of two touching membranes simply add, the partner's share can be subtracted. For
    every measured cell with an interface:

    $$\text{excess} = \frac{\overline{\text{interface}} - \overline{\text{free}} -
    \overline{\text{partner free}}}{\text{SD}(\text{free})}$$

    Here *partner free* is, at every interface point, the mean free-membrane intensity of the partner cell that
    touches there, averaged over the interface. A partner without any free membrane is left out; if no partner of the
    cell has a free membrane, the excess is empty. The control windows have no partner, so their excess is the control
    itself. Files: the column `excess stdDevs above both membranes` of the membrane table, `excess profile` in the
    average profile, and `median excess above both membranes` in `Membrane pairs.tsv`.

    **The check it fails.** For contacts between cells of the same type, where both carry the marker and nothing is
    expected, the excess should be about 0. On real data it was clearly negative for every band placement (outline,
    per type, per cell), every width tried (3, 5 and 9 px) and both ways of reading across the band (*maximum* and
    *sum*): the subtraction takes away more than the partner adds. Centring the band on the membrane even made it more
    negative, because the free membrane is then measured at its peak. So it is not a matter of where the band lies.
    Likely reasons, not proven: the images are 2D projections of 3D cells, so a contact that is tilted is not seen as
    two membranes side by side, and cells that overlap in depth blur differently than free membranes; and a membrane
    may carry less marker per area where it is pressed against another cell.

    **How to read it, if at all.** A positive excess means more marker at the contact than both membranes carry
    elsewhere, whichever cell holds it. Because it tends to be too low, a negative value is not evidence of depletion.
    Check the same-type contacts of your own data in the *All pairs* table first. The comparison of conditions uses
    the log2 ratio to the own membrane.

??? note "Technical: across the band — maximum or sum"

    The dashboard's switch *across the band* chooses how the band becomes one value per point. *maximum* (the
    default, as the macro and the tables) keeps the brightest pixel across the band. *sum* adds all pixels across the
    band. At an interface both membranes lie in the band; only the sum adds them up, which is what *SD above both
    membranes* assumes. The sum also counts the inside of the cell within the band, so it depends more on the band
    width. Everything in the membrane card follows the switch; the tables of the output folder keep the maximum.

## The control

Cells without an interface get **control windows**: *Control windows per cell* (default 10) windows at random places
on their free membrane. Each window is as long as the median interface and is scored like an interface, with both
scores. The histogram of the controls shows how much a stretch of free membrane varies by chance in this dataset; its
log2 ratio lies around 0.

In the dashboard, both histograms share one scale. For the log2 ratio it starts at −3 (one eighth) and reaches at least
2 (four times). For the SD scores it starts at −5 SD, as in the files `Histogram_interface_stdDevs_above_mean -
<marker>.tsv` and `Histogram_control_stdDevs_above_mean - <marker>.tsv` (which stay in SD, as the macro's). If
*Histogram maximum* and *Histogram bin width* are left empty, the scale is automatic. The maximum is the 99th
percentile of the interface values or of the control windows, whichever is higher, but at least 5 SD (2 in log2), and
it is rounded up to a whole number of bins. The bin width is the smallest of 0.1, 0.2, 0.25, 0.5, 1, 2.5, 5, 10, and
so on, that gives at most 60 bins. Values above the maximum are not counted. The two settings apply to the SD scale;
you can change either scale in the dashboard without running the analysis again.

**The control windows are not a control for the interface itself.** Because they lie on free membrane, they cannot
show what touching does to the measurement apart from the marker: two stained membranes side by side, a partner's
signal blurring into the interface, and the outline lying differently at a contact. Therefore:

- **Compare pairs**, not an interface with the control windows. Compare the same measured cells with partners that do
  carry a receptor or ligand and with partners that do not (for example *Immune A+* against *Immune negative*). Or
  compare with contacts between cells of the measured kind (*log2 ratio above same-type contacts*, or set *Interface with* to
  the measured type itself). The dashboard shows all pairs at once (*All pairs* in the membrane card). In the example
  below, the tumor cells' membrane is enriched mainly at contacts with *Subtype A+* cells.

![All pairs of the membrane card](../images/dash_pairs.png)

- A small shift (a ratio below about 1.5, a log2 ratio below about 0.6, or 1–2 SD) is within what the effects above can produce; a large effect that
  is specific to one pair is not.
- Check in the dashboard that the band lies on the membrane (*Cell* card, *membrane band*; membrane card, *Across the
  membrane*).
- A p-value from thousands of cells against thousands of windows says nothing about biological replicates. Compare
  datasets (donors, replicates) instead. A change between conditions is tested across the donors in the overview; see
  [Comparing conditions](comparing-conditions.md).

??? note "Technical: comparing pairs, and why there is no reference channel"

    **The pair table.** For every selection under *Membrane of* (rows) and every cell type or subtype under
    *Interface with* (columns), the dashboard computes the median score of the measured cells that have an interface
    with that partner, and the number of such cells. The median is computed exactly as for the chosen pair, but
    without the control windows. A pair in which the measured cells or the partners do not carry the marker is the
    geometry-matched control for a pair in which they do: the same cells, the same kind of contact and the same
    segmentation. Contacts between cells of the same type show what two stained membranes side by side produce. Each
    entry carries the marker check of its pair: "!" when the partners carry part of the marker or as much as the
    measured cells, "!!" when they carry more, and "–" when the measured cells hardly carry it.

    **No reference channel.** A second channel, measured like the marker of interest but not expected to be enriched,
    seemed a natural control. The cells' own tracker dye was tested as such a reference. Because the dye fills the
    cell, it fills the band at a contact but not at the free edge, and it showed 4.8 SD (at 5 px) and 2.7 SD (at 9 px)
    of "enrichment". A reference must be stained in the same way as the marker (a membrane stain on the same cells),
    which a panel rarely includes. The pair comparison needs no extra channel.

??? note "Technical: the control windows"

    The windows are drawn on the free membrane only, which is read as one closed line, and each window is scored
    against the rest of that free membrane. The macro used one window of 25 px per membrane; many windows of the
    typical interface length give a better estimate of the null distribution. The windows are random and not seeded,
    so the control differs slightly from run to run (the dashboard uses a fixed seed). All windows are listed in
    `Control windows stdDevs above mean - <marker>.tsv`, with both scores.

## Where the marker is expressed

**The marker of interest should be on the measured cells, not on their partners.** At an interface, the membranes of
the two cells lie against each other. If the partner also carries the marker, its membrane is measured together with
that of the measured cell, and the interface looks enriched even when both cells are stained evenly.

The analysis checks this. It compares the median free-membrane intensity of the partners with that of the measured
cells. The result is reported in the log, in `run.json`, in the overview (*MOI partners / measured*) and as a badge in
the dashboard:

| Partners / measured | Badge | Meaning |
|---|---|---|
| below 25 % | ✓ green | The marker is mainly on the measured cells. |
| 25–80 % | ! amber | The partners carry part of the marker, so part of the interface signal may come from them. |
| 80–125 % | ! amber | The marker is on both cell types at a similar level. Two stained membranes side by side look like enrichment; use *log2 ratio above same-type contacts*, or compare with contacts between cells of the same type. |
| 125 % or more | ! red | The marker is brighter on the partners. Check the pair (*Membrane of*, *Interface with*). |
| (any) | – gray | The measured cells hardly carry the marker (below 25 % of the brightest cell type), so their enrichment is not informative. |

??? note "Technical: the marker check, exactly"

    Every cell gets a **free-membrane intensity**: the mean of its membrane profile where its neighbor track shows no
    cell (background-subtracted). Without the membrane analysis, the mean of the whole cell is used instead.
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

If the badge is amber or red, you can choose the pair the other way round (if the marker is on the partners), use
*log2 ratio above same-type contacts*, compare with a pair whose partners do not carry the marker, or treat the result as an
upper limit.

**Look at the subsets, on both sides.** The default pair (all cells of the type that carries the marker, with all
cells of the other type) can hide an effect that only one subset shows. In one test, a receptor that was strongly
polarized on *A+* immune cells gave a median of only 0.2 SD over all immune cells, and a ligand that was enriched only
where *B+* cells touched gave 0.7 SD with all immune cells. The *All pairs* table in the dashboard and
`Membrane pairs.tsv` give every combination of *Membrane of* and *Interface with* at once. Read them before you conclude
that there is no enrichment.

??? note "Technical: inside only, and why it is no cure"

    *Straighten inside the measured cell only* reads the band only where it lies inside the measured cell, so that it
    does not reach into the partner. This is a **sensitivity test, not a separation**. Touching membranes are
    nanometers apart, far below the resolution of the ImageStream, so the blur mixes both membranes in every pixel
    near the interface, whichever pixels are read. If the enrichment drops clearly with *inside only*, part of it may
    come from the partner side. If it stays the same, that is reassuring but not proof.

    *Inside only* can also cut away the measured cell's **own** membrane at a contact. The border between two masks
    need not lie on the membrane, and on real data the measured cells' membrane at a contact lay partly on the
    partner's side of that border. A drop with *inside only* is therefore not by itself proof that the signal came
    from the partner. *Across the membrane* in the dashboard shows where the membrane lies at the contacts.

## Several markers of interest

You can tick more than one channel as *MOI* on the import page (`channels.marker_of_interest` is a list). Every marker
of interest is analyzed in the same way, on the same outlines and neighbor tracks. Each marker gets its own pair: by
default, the cell type that carries that marker against the other type. A marker on the tumor cells and a marker on
the immune cells are therefore each measured on the cells that carry them.

All markers of interest are treated alike, also when there is only one: the name of the marker is part of the names
of its files and columns. The membrane files get " - " and the name before the extension, for example
`Membrane intensity analysis - Marker.tsv`. The per-event columns and the overview columns get the name in brackets,
for example `mean tumor interface (Marker)`. `Membrane pairs.tsv` stays one file, with the marker in its first column.
See [The output folder](../guide/output-folder.md).

In the dashboard, the menu *Marker* in the membrane card switches between the markers. The comparison of conditions
tests each marker separately (see [Comparing conditions](comparing-conditions.md)).

## The per-event table

Independently of the profiles, every event gets the mean of each marker of interest at the interface and on the rest
of the membrane, for each cell type. These values are in `All cell combinations.tsv`, and for the events with
one tumor cell and one immune cell also in `Single Tumor Single Immune cell.tsv`, as in the macro. They are measured
in the membrane band of each cell type (*Membrane / interface width*), not in the straightened profiles.

This is the macro's measure, and it has none of the safeguards described above. It takes all partners of the other
type together, has no control, and does not follow the chosen membrane pair. A partner that carries the marker, or two
stained membranes side by side, raise it as much as a real enrichment does (in a test without any enrichment it was
1.27). The overview column *median … membrane ratio, doublets* shows this measure for the cell type that carries the
marker, and the dashboard shows it in the *Events* table and in the sort *membrane ratio per event*. Read it together with
the marker check, and base your conclusions on the pair comparisons.
