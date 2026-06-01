"""Shared configuration helpers for InsertDiversity scripts."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional


def available_threads() -> int:
    """Return the number of hardware threads available to this process."""
    return os.cpu_count() or 1


def resolve_thread_count(config_value: Optional[Any] = None, workload_size: Optional[int] = None) -> int:
    """Resolve a user-configured thread count.

    By default, InsertDiversity uses every hardware thread reported by the host.
    Users can override that with a positive integer in the consolidated config.
    The workload_size argument is accepted for callers that want to report
    job counts, but the default remains the maximum available hardware threads.
    """
    max_available = available_threads()

    if config_value in (None, "", "auto", "max"):
        requested = max_available
    else:
        requested = int(config_value)
        if requested < 1:
            raise ValueError("threads/max_workers must be a positive integer, 'auto', or omitted")

    return max(1, requested)


def configured_threads(cfg: Dict[str, Any], workload_size: Optional[int] = None) -> int:
    """Read thread settings from a script config.

    Prefer the new 'threads' key, but continue accepting legacy 'max_workers'.
    """
    value = cfg.get("threads", cfg.get("max_workers"))
    return resolve_thread_count(value, workload_size=workload_size)


def load_script_config(config_path: Path, script_name: str) -> Dict[str, Any]:
    """Load either a script-specific config or one step from pipeline_config.json.

    A consolidated config is detected by the presence of a top-level 'scripts'
    object. The named script's nested 'config' block is returned, with top-level
    shared keys such as 'threads' inherited unless the step overrides them.
    """
    with config_path.open("r") as handle:
        cfg = json.load(handle)

    scripts = cfg.get("scripts")
    if not isinstance(scripts, dict):
        return cfg

    if script_name not in scripts:
        raise ValueError(f"Consolidated config does not contain scripts.{script_name}")

    step = scripts[script_name]
    if not isinstance(step, dict) or not isinstance(step.get("config"), dict):
        raise ValueError(f"scripts.{script_name} must contain a config object")

    step_cfg = dict(step["config"])
    for shared_key in ("threads", "max_workers"):
        if shared_key in cfg and shared_key not in step_cfg:
            step_cfg[shared_key] = cfg[shared_key]
    return step_cfg
