"""Command-line interface for SMILES2XYZ."""

import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Final, NoReturn

from smiles2xyz import ConversionError, smiles_to_xyz


def main() -> None:
    """Main entry point for the CLI."""
    parser: ArgumentParser = create_parser()
    args: Namespace = parser.parse_args()

    try:
        xyz_content: str = smiles_to_xyz(
            args.smiles, add_hydrogens=args.hydrogens, optimize=args.optimize
        )

        if args.output:
            output_path: Path = Path(args.output)
            output_path.write_text(xyz_content)
            print(f"XYZ coordinates written to {output_path}", file=sys.stderr)
        else:
            print(xyz_content)

    except ConversionError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(2)


def create_parser() -> ArgumentParser:
    """Create and configure the argument parser.

    Returns:
        Configured ArgumentParser instance.
    """
    parser: ArgumentParser = ArgumentParser(
        prog="smiles2xyz",
        description="Convert SMILES strings to XYZ coordinate format",
    )

    parser.add_argument(
        "smiles",
        type=str,
        help="SMILES string to convert (e.g., 'CCO' for ethanol)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file path (default: print to stdout)",
    )

    parser.add_argument(
        "--no-hydrogens",
        dest="hydrogens",
        action="store_false",
        default=True,
        help="Do not add explicit hydrogen atoms",
    )

    parser.add_argument(
        "--no-optimize",
        dest="optimize",
        action="store_false",
        default=True,
        help="Skip geometry optimization",
    )

    return parser


if __name__ == "__main__":
    main()
