# Settings files and the command line

Everything the import page does can be done from a terminal, with a **project file** (several exports) or a
**settings file** (one dataset). Both are [TOML](https://toml.io) text files.

## Commands

```powershell
uv run imagestream-cli gui                          # the import page (also: uv run imagestream)
uv run imagestream-cli run-project D:\data\project  # all datasets of a project (project.toml or its folder)
uv run imagestream-cli run settings.toml            # one dataset with a settings file
uv run imagestream-cli run settings.toml --max-events 50
uv run imagestream-cli example-settings > settings.toml   # a settings file with every default
```

## A settings file

One dataset: the sections `[input]`, `[channels]`, `[segmentation]`, `[classification]` and `[membrane]`. Relative
paths are relative to the settings file. An unknown setting is an error, so a typing error is not silently ignored.
[All settings](settings.md) lists them.

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
marker_of_interest = 'Marker'
tumor = 'Tumor'
immune = ['Immune']
functional = ['Subtype A', 'Subtype B']
cell_of_interest = 'Tumor cell'

[segmentation]
cellpose_python = 'D:/envs/cellpose3/Scripts/python.exe'
diameter = 20.0

[membrane]
analysis = true
```

The run writes the complete settings it used to `settings.toml` in the dataset's folder, with the channels by number:
a run can be repeated exactly with that file.

## A project file

`project.toml` as the import page writes it: the panel (`[channels]`, by name), the shared settings, and one `[[zip]]`
per export with its datasets:

```toml
[channels]
brightfield = 'BF'
marker_of_interest = 'Marker'
tumor = 'Tumor'
immune = ['Immune']
functional = ['Subtype A', 'Subtype B']
cell_of_interest = 'Tumor cell'

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
```

??? note "Technical: older files"

    Settings files from before November 2026 name the immune marker `t_cell = 4` and the cell of interest
    `'T cell'`; they are read as `immune = [4]` and `'Immune cell'`. Older project files with roles per export, by
    number, are converted to roles by name when all exports agree (see
    [The import page](../guide/import-page.md#2-ideas-exports)).
