# Questions and problems

## Setting up

**The run button stays grey.** The list above it says what is missing: an output folder, a pixel size, or a channel for
a role in one of the exports. A role is a tick box next to a named channel; an export whose channel has another name
(a typing difference) does not get the role.

**"Cellpose failed".** The log shows Cellpose's own message. Check *Python with Cellpose*: it must be the `python.exe`
of an environment in which Cellpose is installed. *Use the GPU* off runs on the processor (much slower).

**An export has objects without all channels.** IDEAS can export objects with missing channels; they are left out,
and the log says how many. *Complete* on the import page counts the objects that have all needed channels.

## Results

**Cells are missed or merged.** Check *Cell diameter (px)* first: Cellpose works best when it matches the cells. For
cells that stain the membrane rather than the cytoplasm, set *Tumor / immune markers* to *membrane*. *Montage columns*
does not change much; leave it at 0.

**Tumor and immune cells are swapped, or mixed.** The k-means needs both cell types and markers that separate them.
Look at *Cell classes* in the dashboard: two clouds should be visible. With very few cells of one type, the
clustering can split the other type instead.

**There are many artefacts.** The automatic thresholds are liberal, so real cells can be positive in all subtype
markers by chance. Look at the histograms in the dashboard; try *MaxEntropy strict*, set the thresholds by hand, or
raise the number of markers in the artefact rule. See [Cell types](../analysis/classification.md#artefacts).

**The partner preference is close to 1 everywhere.** Most contacts may be in events of two cells, which cannot show a
preference and are left out; the dashboard says how many pairs remain. With few events of three or more cells, the
ratios are uncertain. See [Cell contacts](../analysis/interactions.md).

**The membrane badge is amber or red.** The partners carry the marker of interest too, and the interface may show their
membrane. See [Where the marker is expressed](../analysis/membrane.md#where-the-marker-is-expressed).

**The control is not centred on 0.** The control windows score random stretches of free membrane; with uneven
staining along the membrane, the control spreads. The comparison that matters is the interface against the control
of the same dataset.

**The numbers differ slightly between two runs.** The control windows are random (not seeded); everything else is
the same for the same settings. The dashboard's control uses a fixed seed.

## The dashboard

**The dashboard is slow with many events.** Only the visible events are drawn; the plots draw every cell. Filter the
events, or analyse fewer events (*Only the first N events*) to test settings.

**Can I send a dashboard to a colleague?** Yes: send the dataset's folder (`dashboard.html` with `dashboard_files/`).
It opens in any modern browser, without the program and without internet.
