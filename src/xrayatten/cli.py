"""Command line interface for xrayatten."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .local_nist import validate_data


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="xrayatten")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("validate-data", help="Validate local NIST v1.4 manifest and coefficient files")
    run_parser = subparsers.add_parser("run", help="Run a YAML workflow")
    run_parser.add_argument("config", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "validate-data":
        result = validate_data()
        print(
            "validated local NIST data: "
            f"{result['elements']} elements, version {result['version']}, manifest {result['manifest_sha256']}"
        )
        return 0
    if args.command == "run":
        from .workflows import run_config

        result = run_config(args.config)
        print(f"workflow {result.workflow} wrote {len(result.output_paths)} files to {result.output_dir}")
        return 0
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
