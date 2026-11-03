# Questions and problems

## Setting up

**The run button stays gray.** The list above the button shows what is missing: an output folder, a pixel size, or a
channel for one of the roles in one of the exports. A role is a tick box next to a named channel. If a channel in an
export has a different name (for example because of a typing difference), it does not get the role.

**"Cellpose failed".** The log shows Cellpose's own error message. Check *Python with Cellpose*: it must point to the
`python.exe`, or the folder, of an environment in which Cellpose is installed. The chip next to *Segmentation
(Cellpose)* on the import page shows whether Cellpose and a GPU were found there. If *Use the GPU* is off, Cellpose runs on the processor,
which is much slower.

**An export has objects without all channels.** IDEAS can export objects with missing channels. These objects are
left out, and the log says how many there are. On the import page, the *Complete* column counts the objects that have
all the channels that are needed.

## Results

**Cells are missed or merged.** Check *Cell diameter (px)* first: Cellpose works best when it matches the size of the
cells. If the markers stain the membrane rather than the cytoplasm, set *Tumor / immune markers* to *membrane*.
Changing *Montage columns* does not make much difference; leave it at 0.

**Tumor and immune cells are swapped, or mixed.** The k-means clustering needs both cell types to be present, and
markers that separate them. Look at *Cell classes* in the dashboard, where two separate clouds should be visible. If
there are very few cells of one type, the clustering may split the other type into two instead.

**There are many artefacts.** The automatic thresholds are liberal, so real cells can be positive in all subtype
markers by chance. Look at the histograms in the dashboard. Then try *MaxEntropy strict*, set the thresholds by hand,
or, if the artefact rule was set to fewer than all subtype markers, raise that number (0 switches the rule off). See [Cell types](../analysis/classification.md#artefacts).

**A rare subtype appears that should not exist** (for example *A+ C+* in a few cells). The automatic threshold of one
of the markers probably lies within the noise of the negative cells. Look at the histogram of that marker in the
dashboard, drag the slider into the gap between the two peaks, and click *Apply to the project…*. See
[Cell types](../analysis/classification.md#immune-subtypes).

**I only want cell types and contacts.** Leave the *MOI* role unticked and *Membrane profiles and enrichment* off.
The run then does not need a marker of interest.

**I have several markers of interest.** Tick *MOI* for each of them. Each is measured on the cell type that carries it,
and the dashboard's membrane card gets a *Marker* menu. See
[Several markers of interest](../analysis/membrane.md#several-markers-of-interest).

**Does a subtype prefer the tumor cells?** Look at the *among immune cells* tab of the *Interactions* card in the
dashboard. It compares the subtypes with each other within events, with a 95 % interval. See
[Cell contacts](../analysis/interactions.md#subset-preference-which-immune-cells-touch).

**The partner preference is close to 1 everywhere.** Most contacts may be in events of only two cells. Such events
cannot show a preference and are left out; the dashboard says how many pairs remain. If there are few events of three
or more cells, the ratios are uncertain. See [Cell contacts](../analysis/interactions.md).

**The membrane badge is amber or red.** The partner cells carry the marker of interest too, so the interface may show
their membrane rather than that of the measured cells. By default the measured cells are the type that carries the
marker, so a red badge means that *Membrane of* or *Interface with* was set by hand. See
[Where the marker is expressed](../analysis/membrane.md#where-the-marker-is-expressed).

**The control is not centered on 0.** The control windows score random stretches of free membrane. If the staining
is uneven along the membrane, the control values spread out. The comparison that matters is the interface against the
control of the same dataset.

**The numbers differ slightly between two runs.** The control windows are placed at random (without a seed);
everything else is the same when the settings are the same. The control in the dashboard uses a fixed seed.

## The dashboard

**The dashboard is slow with many events.** Only the events that are visible are drawn, but the plots draw every
cell. Filter the events, or analyze fewer events (*Only the first N events*) while you test settings.

**Can I send a dashboard to a colleague?** Yes. Send the folder of the dataset (`dashboard.html` together with
`dashboard_files/`). It opens in any modern browser, without the program and without an internet connection.
