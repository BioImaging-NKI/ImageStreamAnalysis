# ImageStream Analysis

Analysis of cell–cell interactions in imaging flow cytometry (Amnis/Cytek ImageStream) data: segmentation of the
cells in every event image, classification into tumor cells and immune cells with immune subtypes, the cells that
touch, and the enrichment of a marker of interest at the interface between two cell types.

ImageStream Analysis reads IDEAS exports directly, runs the analysis of many exports with one panel of markers, and
writes for every dataset an interactive dashboard and tables, and for the project an overview.

<div class="grid cards" markdown>

- **[Installation](start/install.md)**: uv, the program, and a Python environment with Cellpose.
- **[Your first analysis](start/first-analysis.md)**: a walk-through with the example data, from the IDEAS export to
  the dashboard.
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

ImageStream Analysis for Python is a port of the Fiji macros of
[BioImaging-NKI/ImageStreamAnalysis](https://github.com/BioImaging-NKI/ImageStreamAnalysis): the two-cell-type
analysis with membrane enrichment (v1.8) and the classification of immune subtypes (ClassifyCellTypes 2.8). Where the
port differs from the macros, and why, is listed in [From the Fiji macros](reference/fiji.md).

Authors: Bram van den Broek and Rolf Harkes, The Netherlands Cancer Institute (NKI), Amsterdam. License:
GPL-3.0-or-later.
