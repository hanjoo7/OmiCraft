"""Exclusive GPU leases for independent design jobs on devices 4–7."""

import sys as _sys
if _sys.platform == "win32":
    import types as _types
    _sys.modules.setdefault("fcntl", _types.ModuleType("fcntl"))

import fcntl
import os
import subprocess
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from threading import Event


@dataclass(frozen=True)
class GPUDevice:
    physical: int
    logical: int
    free_mib: int


def available_devices():
    try:
        output = subprocess.run(
            ['nvidia-smi', '--query-gpu=index,uuid,memory.free,utilization.gpu',
             '--format=csv,noheader,nounits'], capture_output=True, text=True, check=True, timeout=10).stdout
        processes = subprocess.run(
            ['nvidia-smi', '--query-compute-apps=gpu_uuid,pid', '--format=csv,noheader,nounits'],
            capture_output=True, text=True, check=True, timeout=10).stdout
        occupied = {line.split(',')[0].strip() for line in processes.splitlines() if ',' in line}
        mask = os.environ.get('CUDA_VISIBLE_DEVICES')
        visible = [v.strip() for v in mask.split(',')] if mask is not None else None
        devices = []
        for line in output.splitlines():
            index, uuid, free, utilization = [part.strip() for part in line.split(',')]
            index, free, utilization = int(index), int(free), int(utilization)
            if index not in {4, 5, 6, 7} or free < 24000 or utilization >= 25 or uuid in occupied:
                continue
            logical = index
            if visible is not None:
                logical = next((i for i, value in enumerate(visible)
                                if value == str(index) or (value.startswith('GPU-') and uuid.startswith(value))), None)
                if logical is None:
                    continue
            devices.append(GPUDevice(index, logical, free))
        return sorted(devices, key=lambda gpu: (-gpu.free_mib, gpu.physical))
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        raise ValueError('GPU 4–7 상태를 확인할 수 없습니다.') from exc


def lock_directory():
    return Path(tempfile.gettempdir()) / ('omicraft-gpu-' + str(os.getuid()))


@contextmanager
def reserve_gpu(*, timeout=1800, poll_seconds=2):
    root = lock_directory()
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    deadline = time.monotonic() + timeout
    while True:
        for device in available_devices():
            handle = (root / f'gpu-{device.physical}.lock').open('a')
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                handle.close()
                continue
            try:
                from .resource_usage import gpu_allocation
                with gpu_allocation(device.physical):
                    yield device
                return
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)
                handle.close()
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise ValueError('GPU 4–7 할당 대기 시간이 초과됐습니다.')
        Event().wait(min(poll_seconds, remaining))
