# Settings files and the command line

Everything that the import page does can also be done from a terminal. You then use either a **project file**, for
several exports, or a **settings file**, for one dataset. Both are [TOML](https://toml.io) text files.

## Commands

```powershell
uv run imagestream-cli gui                          # the import page (also: uv run imagestream)
uv run imagestream-cli shortcut                     # a shortcut to the import page on the Desktop and in the Start menu
uv run imagestream-cli run-project D:\data\project  # all datasets of a project (project.toml or its folder)
uv run imagestream-cli run-project D:\data\project --dataset "Donor 2"   # only this dataset (repeat for several)
uv run imagestream-cli overview D:\data\project     # rewrite the overview tables, without analyzing (after new conditions)
uv run imagestream-cli run settings.toml            # one dataset, with a settings file
uv run imagestream-cli run settings.toml --max-events 50
uv run imagestream-cli example-settings > settings.toml   # a settings file with every default
```

With `--dataset`, `run-project` runs only the named datasets; the others keep their results, and the overview is
written again from all of them. With pools per experiment (`classification.pool = "experiment"`), the other datasets
of the same experiment run too, because they are classified together.

## A settings file

A settings file describes one dataset. It has the sections `[input]`, `[channels]`, `[segmentation]`,
`[classification]` and `[membrane]`. Relative paths are taken relative to the folder of the settings file. An unknown
setting is an error, so that a typing error is not silently ignored. [All settings](settings.md) lists every setting.

```toml
[input]
image = 'D:/data/example_export.zip'
dataset = 'Donor 1'
channels = [1, 2, 3, 4, 7, 8]
channel_names = ['BF', 'Marker', 'Tumor', 'Immune', 'Subtype A', 'Subtype B']
pixel_size = 0.5
output = 'results'

[channels]
brightfield = 'BF'           # a channel by name (or by number: brightfield = 1)
marker_of_interest = ['Marker']  # one or more, or [] for none
tumor = 'Tumor'
immune = ['Immune']
functional = ['Subtype A', 'Subtype B']

[segmentation]
cellpose_python = 'D:/envs/cellpose3/Scripts/python.exe'
diameter = 20.0

[membrane]
analysis = true
```

Each run writes the complete settings that it used to `settings.toml` in the folder of the dataset, with the channels
given by number. You can use that file to repeat the run exactly.

## A project file

The import page writes the project file `project.toml`. It holds the panel (the `[channels]` section, with the
channels given by name), the settings that all datasets share, and one `[[zip]]` section for each export, with its
datasets:

```toml
[channels]
brightfield = 'BF'
marker_of_interest = ['Marker']
tumor = 'Tumor'
immune = ['Immune']
functional = ['Subtype A', 'Subtype B']

[classification]
pool = 'experiment'

[[zip]]
path = 'D:/data/example_export.zip'
pixel_size = 0.5
channels = [1, 2, 3, 4, 7, 8]
channel_names = ['BF', 'Marker', 'Tumor', 'Immune', 'Subtype A', 'Subtype B']

[[zip.dataset]]
name = 'Donor 1'
experiment = 'Example'
include = true
condition = 'control'      # optional, to compare conditions: the condition
replicate = 'Donor 1'      # and the independent unit (donor) it comes from
```

The reference condition, with which the other conditions are compared, goes in a section of its own. If it is empty
or absent, the first condition in the file is used:

```toml
[comparison]
reference = 'control'
```

??? note "Technical: older files"

    Settings files from before November 2026 name the immune marker `t_cell = 4`; it is read as `immune = [4]`. They
    give one marker of interest (`marker_of_interest = 2` or `'Marker'`, or `''` for none), which is read as a list
    (`[2]`, `['Marker']` or `[]`). Their `cell_of_interest` (*Marker of interest expressed in*) is ignored with a
    note in the log: the measured cells now follow from the data (see
    [A pair of cell types](../analysis/membrane.md#a-pair-of-cell-types)). Older project files, which set the roles
    for each export by channel number, are converted to roles by name when all exports agree (see
    [The import page](../guide/import-page.md#2-ideas-exports)).
