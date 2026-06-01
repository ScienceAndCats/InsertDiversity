#!/usr/bin/env python3
"""Run InsertDiversity scripts from one consolidated JSON config."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

from config_helpers import run_cli


COMMANDS = {
    "bun_extract": [sys.executable, "bun_extract.py", "--config"],
    "barcode_from_bun_csv": [sys.executable, "barcode_from_bun_csv.py", "--config"],
    "ordered_barcode_sample_analysis": [
        sys.executable,
        "ordered_barcode_sample_analysis.py",
    ],
}


def load_pipeline_config(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Pipeline config not found: {path}")

    with path.open("r") as handle:
        cfg = json.load(handle)

    if "scripts" not in cfg or not isinstance(cfg["scripts"], dict):
        raise ValueError("Pipeline config must contain a 'scripts' object.")

    return cfg


def selected_script_names(cfg: Dict[str, Any]) -> List[str]:
    scripts = cfg["scripts"]
    run_order = cfg.get("run_order") or list(scripts.keys())

    missing = [name for name in run_order if name not in scripts]
    if missing:
        raise ValueError("run_order includes scripts not present in scripts: " + ", ".join(missing))

    return [name for name in run_order if scripts[name].get("enabled", False)]


def command_for(name: str, step: Dict[str, Any], config_path: Path) -> List[str]:
    command = step.get("command")
    if command is None:
        command = COMMANDS.get(name)

    if not command:
        script = step.get("script")
        if not script:
            raise ValueError(
                f"Script '{name}' needs either a known command, a command array, or a script path."
            )
        script_path = Path(script)
        if not script_path.is_absolute():
            script_path = Path(__file__).resolve().parent / script_path
        if not script_path.is_file():
            raise FileNotFoundError(f"Script file for '{name}' not found: {script_path}")
        command = [sys.executable, script, "--config"]

    if not isinstance(command, list) or not all(isinstance(x, str) for x in command):
        raise ValueError(f"Script '{name}' command must be a list of strings.")

    return [*command, str(config_path)]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run enabled InsertDiversity analysis steps from a consolidated JSON config."
    )
    parser.add_argument(
        "--config",
        default="pipeline_config.json",
        help="Path to consolidated pipeline JSON config (default: pipeline_config.json).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the enabled commands without running them.",
    )
    args = parser.parse_args()

    pipeline_config_path = Path(args.config)
    if not pipeline_config_path.is_absolute():
        pipeline_config_path = pipeline_config_path.resolve()
    cfg = load_pipeline_config(pipeline_config_path)
    enabled_names = selected_script_names(cfg)

    if not enabled_names:
        print("No scripts are enabled. Set at least one scripts.<name>.enabled value to true.")
        return

    repo_root = Path(__file__).resolve().parent
    planned_commands = []

    for name in enabled_names:
        step = cfg["scripts"][name]
        planned_commands.append((name, command_for(name, step, pipeline_config_path)))

    if args.dry_run:
        print("Enabled pipeline steps:")
        for name, command in planned_commands:
            print(f"  {name}: {' '.join(command)}")
        return

    for name, command in planned_commands:
        print(f"\n=== Running {name} ===")
        print("Command:", " ".join(command))
        try:
            subprocess.run(command, cwd=repo_root, check=True)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Could not start pipeline step '{name}'. Missing executable or file: {exc.filename}"
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"Pipeline step '{name}' failed with exit code {exc.returncode}: {' '.join(command)}"
            ) from exc

    print("\nPipeline complete.")


if __name__ == "__main__":
    run_cli(main)
