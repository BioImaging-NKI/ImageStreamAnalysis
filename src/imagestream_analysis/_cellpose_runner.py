"""Run Cellpose (3 or 4) on one 2-D image. Executed as a standalone script, possibly by the Python of another
environment, so it must not import imagestream_analysis.

Usage: python _cellpose_runner.py input.npy output.npy '{"model": "cyto3", ...}'

input.npy: (Y, X) or (Y, X, C) float32. With Cellpose 3, params["channels"] selects the channels as in the Cellpose
CLI (--chan, --chan2; 0 = grayscale, 1-based otherwise); Cellpose 4 (cpsam) uses up to 3 channels as they are.
"""

import json
import sys
import time
from typing import Any

import numpy as np


def run(img: np.ndarray, params: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    """Segment img; returns the masks and info on the run (Cellpose version, device, timing)."""
    from importlib.metadata import version

    from cellpose import models

    major = int(version("cellpose").split(".")[0])
    model_name = params["model"]
    gpu = params.get("use_gpu", True)
    diameter = params.get("diameter") or None
    is_path = any(ch in model_name for ch in "/\\.")
    kwargs = {"flow_threshold": params["flow_threshold"], "cellprob_threshold": params["cellprob_threshold"]}
    channel_axis = 2 if img.ndim == 3 else None
    t0 = time.perf_counter()

    if major >= 4:
        if model_name != "cpsam" and not is_path:
            raise ValueError(f"Cellpose {version('cellpose')} only has the model 'cpsam' (got '{model_name}')")
        model = models.CellposeModel(gpu=gpu, pretrained_model=model_name)
        t1 = time.perf_counter()
        masks = model.eval(img, diameter=diameter, channel_axis=channel_axis, **kwargs)[0]
    else:
        if model_name == "cpsam":
            raise ValueError(f"The model 'cpsam' needs Cellpose 4 (this environment has {version('cellpose')})")
        channels = params.get("channels") or [0, 0]
        if is_path:
            model = models.CellposeModel(gpu=gpu, pretrained_model=model_name)
        elif diameter is None:
            model = models.Cellpose(gpu=gpu, model_type=model_name)  # includes the size model
        else:
            model = models.CellposeModel(gpu=gpu, model_type=model_name)
        t1 = time.perf_counter()
        masks = model.eval(img, diameter=diameter, channels=channels, channel_axis=channel_axis, **kwargs)[0]
    t2 = time.perf_counter()
    device = getattr(model, "device", None) or getattr(getattr(model, "cp", None), "device", None)
    info = {"version": version("cellpose"), "device": str(device), "gpu_requested": bool(gpu),
            "load_s": round(t1 - t0, 2), "eval_s": round(t2 - t1, 2)}
    try:
        import torch

        info["torch"] = torch.__version__
        if torch.cuda.is_available():
            info["gpu_name"] = torch.cuda.get_device_name(0)
    except Exception:  # noqa: BLE001, S110 - info only
        pass
    return np.asarray(masks).astype(np.int32), info

if __name__ == "__main__":
    image = np.load(sys.argv[1])
    masks, info = run(image, json.loads(sys.argv[3]))
    np.save(sys.argv[2], masks)
    print("ISA_INFO " + json.dumps(info))
