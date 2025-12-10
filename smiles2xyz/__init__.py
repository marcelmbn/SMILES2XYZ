"""SMILES2XYZ - Convert SMILES strings to XYZ format.

This package provides functionality to convert SMILES (Simplified Molecular Input
Line Entry System) strings to XYZ coordinate files with explicit hydrogen atoms.
"""

from smiles2xyz.converter import ConversionError, smiles_to_xyz

__version__: str = "0.1.0"
__all__: list[str] = ["smiles_to_xyz", "ConversionError"]
