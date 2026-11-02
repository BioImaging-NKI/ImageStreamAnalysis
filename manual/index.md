# ImageStream Analysis

Which immune cells seek out the tumor cells, and what gathers where they meet? An imaging flow cytometer
(Amnis/Cytek ImageStream) photographs thousands of these encounters per sample. ImageStream Analysis outlines every
cell in every event image, sorts the cells into tumor cells and immune subtypes, finds the ones that touch, and
measures how strongly a marker of interest concentrates at the interface between two cell types.

It takes the IDEAS exports as they are, analyses many of them with one marker panel, and gives every dataset an
interactive dashboard (thresholds and cell-type pair adjustable live) and the project an overview that compares
conditions across donors.

[Try the live demo](https://bioimaging-nki.github.io/ImageStreamAnalysis/demo/dashboard.html){ .md-button .md-button--primary target="_blank" }
[Your first analysis](start/first-analysis.md){ .md-button }

<div class="isa-flow">
<a class="isa-step" href="analysis/input/"><svg viewBox="0 0 96 72" role="img" aria-label="Read the exports"><rect x="30" y="4" width="40" height="34" rx="3" fill="#3a3a3a"/><rect x="22" y="12" width="40" height="34" rx="3" fill="#3a3a3a"/><rect x="14" y="20" width="40" height="34" rx="3" fill="#3a3a3a"/>
  <circle cx="30" cy="37" r="9" fill="#eb6834"/><circle cx="43" cy="33" r="5.5" fill="#2a78d6"/>
  <g stroke="currentColor" stroke-width="1.5" fill="none" opacity=".75"><path d="M73 30h8m-4-4 4 4-4 4"/></g>
  <rect x="85" y="14" width="9" height="9" rx="2" fill="#9a9a9a"/><rect x="85" y="27" width="9" height="9" rx="2" fill="#eb6834"/><rect x="85" y="40" width="9" height="9" rx="2" fill="#2a78d6"/></svg>
<strong>Read the exports</strong><span>IDEAS zips as exported, no conversion step: every object, every channel, named by its marker.</span></a>
<a class="isa-step" href="analysis/segmentation/"><svg viewBox="0 0 96 72" role="img" aria-label="Segment"><rect x="8" y="6" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="19" cy="21" r="6" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="27" cy="17" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><rect x="36" y="6" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="47" cy="21" r="8" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="55" cy="17" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><rect x="64" y="6" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="75" cy="21" r="6" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="83" cy="17" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><rect x="8" y="36" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="19" cy="51" r="8" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="27" cy="47" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><rect x="36" y="36" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="47" cy="51" r="6" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="55" cy="47" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><rect x="64" y="36" width="26" height="28" rx="3" fill="#3a3a3a"/><circle cx="75" cy="51" r="8" fill="none" stroke="#f2f2f2" stroke-width="1.6"/><circle cx="83" cy="47" r="4" fill="none" stroke="#f2f2f2" stroke-width="1.6"/></svg>
<strong>Segment</strong><span>Cellpose outlines the cells of all events at once, on one montage; each cell is measured in every channel.</span></a>
<a class="isa-step" href="analysis/classification/"><svg viewBox="0 0 96 72" role="img" aria-label="Classify"><circle cx="26" cy="34" r="16" fill="#eb6834"/>
  <circle cx="54" cy="18" r="7" fill="#2a78d6"/><circle cx="72" cy="26" r="7" fill="#eda100"/><circle cx="58" cy="44" r="7" fill="#e87ba4"/>
  <circle cx="80" cy="50" r="5" fill="#9a9a9a"/><path d="M76 46l8 8m0-8-8 8" stroke="#f2f2f2" stroke-width="1.4"/></svg>
<strong>Classify</strong><span>Tumor or immune cell by k-means, immune subtypes by a threshold per marker, and doubtful detections set aside as artefacts.</span></a>
<a class="isa-step" href="analysis/interactions/"><svg viewBox="0 0 96 72" role="img" aria-label="Contacts"><circle cx="64.9" cy="27.6" r="8" fill="#2a78d6"/><circle cx="28.1" cy="20.1" r="8" fill="#eda100"/><circle cx="80" cy="56" r="7" fill="#e87ba4"/>
  <circle cx="44" cy="36" r="16" fill="#eb6834"/>
  <path d="M56.8 26.4A16 16 0 0 1 59.9 34.1" stroke="#1baf7a" stroke-width="3.5" fill="none" stroke-linecap="round"/>
  <path d="M30.1 28.0A16 16 0 0 1 36.0 22.1" stroke="#1baf7a" stroke-width="3.5" fill="none" stroke-linecap="round"/>
  <path d="M68 35 78 49" stroke="currentColor" stroke-dasharray="2 3" stroke-width="1.3" opacity=".6"/></svg>
<strong>Contacts</strong><span>Who touches whom, and more often than chance? Each event is compared with its own cells shuffled.</span></a>
<a class="isa-step" href="analysis/membrane/"><svg viewBox="0 0 96 72" role="img" aria-label="Membrane enrichment"><circle cx="46.9" cy="21.4" r="8" fill="#2a78d6"/>
  <circle cx="26" cy="36" r="17" fill="none" stroke="#eb6834" stroke-width="4"/>
  <path d="M36.5 22.6A17 17 0 0 1 42.2 30.7" stroke="#1baf7a" stroke-width="5" fill="none" stroke-linecap="round"/>
  <path d="M60 58H92" stroke="currentColor" stroke-width="1" opacity=".5"/>
  <polyline points="60,55 66,54 70,55 74,40 77,22 80,40 84,54 88,55 92,54" fill="none" stroke="#1baf7a" stroke-width="2"/></svg>
<strong>Membrane enrichment</strong><span>Does the marker gather where the cells meet? The interface against the rest of the membrane and random control windows.</span></a>
<a class="isa-step" href="guide/overview/"><svg viewBox="0 0 96 72" role="img" aria-label="Results"><rect x="6" y="6" width="84" height="60" rx="5" fill="none" stroke="currentColor" stroke-width="1.5" opacity=".8"/>
  <path d="M6 16H90" stroke="currentColor" stroke-width="1" opacity=".5"/>
  <rect x="11" y="21" width="14" height="13" rx="3" fill="#3a3a3a"/><rect x="28" y="21" width="14" height="13" rx="3" fill="#3a3a3a"/><rect x="11" y="37" width="14" height="13" rx="3" fill="#3a3a3a"/><rect x="28" y="37" width="14" height="13" rx="3" fill="#3a3a3a"/>
  <circle cx="18" cy="28" r="3.5" fill="#eb6834"/><circle cx="36" cy="27" r="3" fill="#2a78d6"/><circle cx="17" cy="44" r="3" fill="#eda100"/><circle cx="35" cy="44" r="3.5" fill="#eb6834"/>
  <g fill="#2a78d6" opacity=".9"><rect x="50" y="40" width="5" height="16"/><rect x="57" y="30" width="5" height="26"/><rect x="64" y="35" width="5" height="21"/><rect x="71" y="46" width="5" height="10"/><rect x="78" y="50" width="5" height="6"/></g>
  <path d="M48 24h36M48 30h28" stroke="currentColor" stroke-width="1.5" opacity=".5"/></svg>
<strong>Results</strong><span>A live dashboard per dataset, a project overview, and conditions tested against each other across donors.</span></a>
</div>

<div class="grid cards" markdown>

- **[Installation](start/install.md)**: uv, the program, and a Python environment with Cellpose.
- **[Your first analysis](start/first-analysis.md)**: a walk-through with the example data, from the IDEAS export to
  the dashboard.
- **[Try it in your browser](start/demo.md)**: the dashboard of an example dataset, live, without installing
  anything.
- **[How the analysis works](analysis/input.md)**: every step, with the reasoning behind it in foldable *Technical*
  sections.
- **[All settings](reference/settings.md)**: every setting, its default and its meaning.

</div>

## What it does

1. **Reads the IDEAS export**: one folder per dataset, one image per object and channel. The channels get a name and a
   role: brightfield, tumor marker, immune-cell marker(s), immune subtype markers, and the marker of interest.
2. **Segments the cells** with Cellpose, on a montage of all events, and measures every cell in every channel.
3. **Classifies the cells**: tumor or immune cell (k-means), the immune cells into subtypes (a threshold per marker),
   and detections that are positive in too many subtype markers as artefacts. The datasets of one experiment are
   classified together.
4. **Finds the cells that touch** and compares the contacts with chance: which cell types prefer each other within
   the events.
5. **Measures the marker of interest along the membrane** of one cell type, at its interface with another, against
   the rest of its membrane and against random control windows; and checks that the marker is on the measured cells
   rather than on their partners.
6. **Shows the results** in a dashboard per dataset, in which the subtype thresholds and the pair of cell types can be
   changed live, and in a project overview.

## Background

ImageStream Analysis for Python is a port of the Fiji macros that this repository started with (kept in
[`legacy/`](https://github.com/BioImaging-NKI/ImageStreamAnalysis/tree/main/legacy)): the two-cell-type
analysis with membrane enrichment (v1.8) and the classification of immune subtypes (ClassifyCellTypes 2.8). Where the
port differs from the macros, and why, is listed in [From the Fiji macros](reference/fiji.md).

Authors: Bram van den Broek and Rolf Harkes, The Netherlands Cancer Institute (NKI), Amsterdam. License:
GPL-3.0-or-later.
