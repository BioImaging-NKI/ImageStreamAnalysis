# Maintaining the user manual of ImageStream Analysis

Notes for whoever (person or Claude, with the `update-manual` skill) updates the manual after a change to the program.
The manual: `manual/` (MkDocs Material, `mkdocs.yml`), to be published at
<https://bioimaging-nki.github.io/ImageStreamAnalysis/> by `.github/workflows/manual.yml`. Rule: every noticeable
change to the import page, the dashboard, the analysis or the outputs is checked against the manual in the same piece
of work.

**The manual is public.** It shows only the synthetic example (below), never research data: no sample, donor, marker
or channel names of real data, no numbers measured on them, and no folders of this computer (the screenshot script
replaces the paths shown on the import page).

## Which code goes with which pages

| Changed code | Check these pages |
|---|---|
| `templates/gui.html`, `gui.py` (steps, labels, checks, buttons) | `guide/import-page.md`, `start/first-analysis.md`, `help/faq.md` |
| `config.py` (settings, defaults, their comments), FIELDS in `gui.html` | regenerate `reference/settings.md`; pages that mention the setting (`grep -rn` its name and label) |
| `templates/dashboard.html`, `dashboard.py` | `guide/dashboard.md`, `start/first-analysis.md` (step 7), `start/demo.md` |
| `project.py` (overview, pools, `interactions.tsv`) | `guide/overview.md`, `guide/import-page.md`, `analysis/classification.md` (pools) |
| `ideas.py`, `io.py`, `measure.py` | `analysis/input.md`, `analysis/segmentation.md` |
| `segmentation.py`, `_cellpose_runner.py` | `analysis/segmentation.md`, `start/install.md` |
| `classification.py`, `imagej.py` (thresholds) | `analysis/classification.md`, `guide/dashboard.md` (subtype markers) |
| `interactions.py` | `analysis/interactions.md`, `guide/dashboard.md` (interactions), `guide/overview.md` |
| `membrane.py`, `interfaces.py` | `analysis/membrane.md`, `guide/dashboard.md` (membrane enrichment) |
| `pipeline.py` output files and columns | `guide/output-folder.md` |
| `cli.py`, settings and project files | `reference/command-line.md` |
| a deliberate difference from the macros (SPEC.md) | `reference/fiji.md` |

Every analysis topic has foldable `??? note "Technical: …"` sections for the reasoning and the pitfalls (for example:
the comparison with chance only in events of three or more cells; what the selection in IDEAS means; where the marker
of interest is expressed). New intricacies go there, not into the main text.

Generated, never edit by hand: `src/imagestream_analysis/templates/help.json`, the context help of the import page and
the dashboard (`python tools/manual/make_help.py`, after every change to the guide, analysis, reference or help pages;
`tests/test_help.py` fails when it is stale or a ? button (`data-help` in the templates) points at a heading that
is gone: rename a heading together with its targets). `manual/reference/settings.md` (`python tools/manual/gen_reference.py`, from the
comments of the settings in `config.py` and the labels in FIELDS of `gui.html`). `tests/test_manual.py` fails when it
is stale, when a setting has no comment, or when the site does not build in strict mode.

## The example

- Synthetic data, made by `python tools/manual/make_example_data.py` (deterministic): `data/manual/example_export.zip`
  (two datasets, *Donor 1* and *Donor 2*; channels BF, Marker, Tumor, Immune, Subtype A, Subtype B and an unused one)
  and `truth.csv` (the cells that were drawn). `data/` is not in git.
- The project: `python tools/manual/make_example_project.py` sets it up as on the import page and runs it
  (`data/manual/project/`, about a minute with a GPU). Cellpose: `cellpose_python` in `data/manual/local.toml`
  (a path of this computer, not in git) or the environment variable `ISA_CELLPOSE_PYTHON`.
- Settings of the example (the walk-through quotes them): marker of interest expressed in tumor cells, cell diameter
  20 px, membrane analysis on, pools per experiment ("Example").
- The example zip is also an asset of the `example-data` release, for readers.

## Screenshots

`python tools/manual/screenshots.py [names]` (Playwright with the installed Edge, 1440×900 at scale 1.5). The import
page runs in that process with a state folder of its own; the user's settings are not touched.

| Screenshot | Shows | Pages |
|---|---|---|
| `gui_exports`, `gui_channels` | step 2 of the import page | first-analysis, import-page |
| `gui_settings` | step 3 | first-analysis |
| `gui_overview` | step 5, with the partner-preference matrix | first-analysis, overview |
| `dash_overview` | the top of the dashboard | first-analysis, dashboard |
| `dash_markers`, `dash_classes`, `dash_cell` | right-hand cards | dashboard, classification |
| `dash_help` | the context help panel, opened from the ? of *Subtype markers* | dashboard |
| `dash_interactions`, `dash_contact_filter` | interactions; the gallery filtered on a contact | dashboard, interactions |
| `dash_membrane` | membrane enrichment | dashboard, membrane |

After a change to a page of the program: remake the example project (its dashboards are written by the program), then
the screenshots, and look at every changed image.

## Waiting for the next refresh

Small changes to the program whose screenshots are not retaken yet (the user prefers to retake them together, with
a larger change). Retake them at the next refresh and empty this list.

- 2026-11-01: the import page's analysis settings in rows of two (three where they belong together), with
  *Reuse the previous Cellpose result* and *Thresholds set by hand*: `gui_settings`.
- 2026-11-01: dashboard buttons *Apply to the project…* and *Copy as settings* (was *Copy the thresholds as
  settings*), *Download .tsv* in the tables, *Apply* in the membrane card: `dash_markers`, `dash_membrane`,
  `dash_overview` (if the tables show).

## Refresh, check, publish

1. `python tools/manual/make_example_project.py` (if the analysis or the dashboard changed)
2. `python tools/manual/screenshots.py`
3. `python tools/manual/gen_reference.py` and `python tools/manual/make_help.py`
4. `python tools/manual/make_demo.py` (the live demo: `data/manual/manual-demo.zip`, and `manual/demo/` locally)
5. `python -m pytest tests/test_manual.py` and `mkdocs serve` to look at the pages
6. Publishing (with the user's OK): upload the demo first, then push:
   `gh release upload example-data data/manual/manual-demo.zip data/manual/example_export.zip --clobber`, then
   `git push`. The workflow builds the site with the pinned versions and deploys it; check the live page afterwards
   (fetch it and search for the new images).

GitHub Pages has to be enabled once by an admin of the repository: Settings → Pages → Source: *GitHub Actions*.
