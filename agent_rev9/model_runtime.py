"""Shared subprocess isolation; gpu_device indexes the parent's visible devices."""

import os


def gpu_environment(gpu_device: int, *, use_gpu=True, environment=None):
    if type(gpu_device) is not int or gpu_device < 0:
        raise ValueError("gpu_device must be a nonnegative integer")
    overrides = environment or {}
    if "CUDA_VISIBLE_DEVICES" in overrides:
        raise ValueError(
            "Select GPUs through screening_options.gpu_device, not environment overrides"
        )
    env = os.environ.copy()
    env.update(overrides)
    if not use_gpu:
        env["CUDA_VISIBLE_DEVICES"] = ""
    elif "CUDA_VISIBLE_DEVICES" in env:
        visible = [v.strip() for v in env["CUDA_VISIBLE_DEVICES"].split(",") if v.strip()]
        if gpu_device >= len(visible) or visible[gpu_device] == "-1":
            raise ValueError("Configured GPU index is not visible")
        env["CUDA_VISIBLE_DEVICES"] = visible[gpu_device]
    else:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_device)
    return env
