# Installation

ImageStream Analysis is a Python program. It runs on Windows (the import page uses Windows file dialogs); the
command line also runs on macOS and Linux. Cellpose, the segmentation, preferably runs on a GPU.

## 1. uv

The program and its packages are installed with [uv](https://docs.astral.sh/uv/), a Python package manager. Install
it once, in PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## 2. The program

Download the repository (on GitHub: *Code* → *Download ZIP*, or `git clone`), open a terminal in its folder, and
install:

```powershell
uv sync
```

This makes a Python environment in the folder `.venv`, with the program and the packages it needs.

## 3. Cellpose

The cells are segmented with [Cellpose](https://www.cellpose.org/), version 3 (model `cyto3`, the default) or 4
(model `cpsam`). There are two ways:

=== "A separate environment (recommended)"

    Install Cellpose in an environment of its own, with the PyTorch version for your GPU, following the Cellpose
    instructions. ImageStream Analysis then calls that environment's Python for the segmentation: on the import page
    under *Python with Cellpose*, or in a settings file as `segmentation.cellpose_python`. The Cellpose environment
    can be shared with other programs, and it can be updated without touching this one.

=== "In this environment"

    ```powershell
    uv sync --extra cellpose3
    ```

    or `--extra cellpose4`. Leave *Python with Cellpose* empty. This installs the CPU or CUDA version of PyTorch that
    pip chooses; for a GPU you may have to install the matching PyTorch yourself.

## 4. Start

```powershell
uv run imagestream
```

opens the [import page](../guide/import-page.md) in your browser. `uv run imagestream-cli` is the
[command line](../reference/command-line.md).

??? note "Technical: why Cellpose runs in another process"

    The segmentation runs as a separate process with the Python you name (`_cellpose_runner.py`): the montage of all
    events is handed over as a file, and the labels come back the same way. So Cellpose 3 and 4, with their own
    PyTorch and CUDA versions, never have to fit in the same environment as this program, and a crash of the GPU
    driver ends that process, not the analysis. The run's `run.json` records the Cellpose version, the device (GPU or
    CPU) and the time it took.
