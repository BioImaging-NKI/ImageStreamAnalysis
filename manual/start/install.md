# Installation

ImageStream Analysis is a Python program. It runs on Windows, because the import page uses Windows file dialogs; the
command line also runs on macOS and Linux. Cellpose, which does the segmentation, preferably runs on a GPU.

## 1. uv

The program and the packages it needs are installed with [uv](https://docs.astral.sh/uv/), a Python package
manager. You install uv once, from PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## 2. The program

Download the repository (on GitHub, click *Code* → *Download ZIP*, or use `git clone`). Then open a terminal in its
folder and install the program:

```powershell
uv sync
```

This creates a Python environment in the folder `.venv` that contains the program and the packages it needs.

## 3. Cellpose

The cells are segmented with [Cellpose](https://www.cellpose.org/), either version 3 (with the model `cyto3`, the
default) or version 4 (with the model `cpsam`). There are two ways to install it:

=== "A separate environment (recommended)"

    Install Cellpose in an environment of its own, with the PyTorch version that suits your GPU, following the
    Cellpose instructions. ImageStream Analysis then uses the Python of that environment for the segmentation. You
    set it on the import page under *Python with Cellpose*, or in a settings file as `segmentation.cellpose_python`.
    Other programs can share the same Cellpose environment, and you can update it without touching this program.

=== "In this environment"

    ```powershell
    uv sync --extra cellpose3
    ```

    or `--extra cellpose4`. Leave *Python with Cellpose* empty. This installs whichever version of PyTorch pip
    chooses, for the CPU or for CUDA. To use a GPU, you may have to install the matching PyTorch version yourself.

## 4. Start

```powershell
uv run imagestream
```

This command opens the [import page](../guide/import-page.md) in your browser. `uv run imagestream-cli` starts the
[command line](../reference/command-line.md) instead.

??? note "Technical: why Cellpose runs in another process"

    The segmentation runs as a separate process (`_cellpose_runner.py`), with the Python that you name. The montage
    of all events is handed over to it as a file, and the labels come back the same way. As a result, Cellpose 3 and
    4, each with their own PyTorch and CUDA versions, never have to fit into the same environment as this program,
    and a crash of the GPU driver ends only that process, not the analysis. The `run.json` file of each run records
    the Cellpose version, the device (GPU or CPU) and the time the segmentation took.
