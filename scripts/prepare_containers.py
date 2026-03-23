#!/usr/bin/env python3
"""
prepare_containers.py — Pull Apptainer SIF images that are absent locally.

Reads config/containers.yaml and runs `apptainer pull` for any SIF that does
not yet exist. Safe to run repeatedly (idempotent).

Usage:
    python scripts/prepare_containers.py [containers_yaml]
    python scripts/prepare_containers.py config/containers.yaml --pull
"""

import argparse
import os
import subprocess
import sys

import yaml


def prepare_containers(containers_yaml: str, pull: bool = True) -> None:
    with open(containers_yaml) as fh:
        containers = yaml.safe_load(fh)

    missing = []
    for tool, spec in containers.items():
        sif = spec["sif"]
        if not os.path.exists(sif):
            missing.append((tool, spec))

    if not missing:
        print("[prepare_containers] All SIF files present.")
        return

    for tool, spec in missing:
        sif = spec["sif"]
        uri = spec["uri"]
        if pull:
            print(f"[prepare_containers] Pulling {tool} ({spec['version']}) → {sif}")
            os.makedirs(os.path.dirname(sif), exist_ok=True)
            subprocess.run(["apptainer", "pull", sif, uri], check=True)
        else:
            print(
                f"[prepare_containers] WARNING: SIF not found: {sif}\n"
                f"  Run: apptainer pull {sif} {uri}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pull missing Apptainer SIF images for the BSMN pipeline."
    )
    parser.add_argument(
        "containers_yaml",
        nargs="?",
        default="config/containers.yaml",
        help="Path to containers.yaml (default: config/containers.yaml)",
    )
    parser.add_argument(
        "--pull",
        action="store_true",
        default=True,
        help="Pull missing SIF files (default: True)",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check for missing SIFs; do not pull",
    )
    args = parser.parse_args()

    pull = not args.check_only
    prepare_containers(args.containers_yaml, pull=pull)


if __name__ == "__main__":
    main()
