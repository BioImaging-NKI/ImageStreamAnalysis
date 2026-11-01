# Glossary

Artefact
:   A detection positive in all (or at least N) subtype markers: taken not to be a cell. Left out of the clustering and
    the counts of cells.

Background
:   Per channel, the median over the objects of each object's median intensity. Subtracted from the cell means.

Dataset
:   A folder in an IDEAS export: the objects of one sample. Analysed into a folder of its own.

Event
:   One object of the ImageStream: one image with one or more cells.

Experiment
:   A label that groups datasets: classified together (a pool) and summed in the overview.

Free membrane
:   The part of a cell's membrane that touches no other cell: the baseline of the enrichment.

Interface
:   The part of a cell's membrane within the interface width (*Membrane / interface width*) of a partner cell.

Marker of interest (MOI)
:   The marker whose enrichment at the interface is measured.

Neighbour track
:   For every position along a cell's membrane, the cell that touches it there (or none).

Partner preference
:   Contacts between two cell types observed ÷ expected by chance, within events of three or more cells.

Pool
:   The datasets that are classified together: one clustering and one set of thresholds.

Project
:   Exports, channel names and roles, and settings, saved as `project.toml`, with the results of all datasets.

Subtype marker
:   A marker (role *Subtype*) that is positive or negative per immune cell; the positive markers name the subtype.

Touching
:   Two cells of one event within the interface width of each other.
