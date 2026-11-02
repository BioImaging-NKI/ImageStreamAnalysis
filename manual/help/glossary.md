# Glossary

Artefact
:   A detection that is positive in all subtype markers (or in at least N of them) and is therefore taken not to be a
    cell. Artefacts are left out of the clustering and out of the cell counts.

Background
:   For each channel, the median, over all objects, of the median intensity of each object. The background is
    subtracted from the mean intensities of the cells.

Dataset
:   A folder in an IDEAS export, which holds the objects of one sample. Each dataset is analyzed into a folder of its
    own.

Event
:   One object recorded by the ImageStream: one image that contains one or more cells.

Experiment
:   A label that groups datasets. Datasets with the same experiment label are classified together (as one pool), and
    their contacts are pooled in the overview's partner-preference matrix.

Free membrane
:   The part of a cell's membrane that touches no other cell. It is the baseline for the enrichment.

Interface
:   The part of a cell's membrane that lies within the interface width (*Membrane / interface width*) of a partner
    cell.

Marker of interest (MOI)
:   The marker whose enrichment at the interface is measured.

Neighbor track
:   For every position along a cell's membrane, the cell that touches it at that position (or none).

Partner preference
:   The number of contacts between two cell types that is observed, divided by the number expected by chance, in
    events of three or more cells.

Pool
:   A group of datasets that are classified together, with one clustering and one set of thresholds.

Project
:   The exports, the channel names and roles, and the settings, saved as `project.toml`, together with the results of
    all datasets.

Subtype marker
:   A marker (with the role *Subtype*) that is scored as positive or negative in each immune cell. The positive
    markers give the subtype its name.

Touching
:   Two cells in the same event are touching when they lie within the interface width of each other.
