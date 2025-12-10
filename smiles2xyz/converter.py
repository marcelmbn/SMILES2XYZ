"""SMILES to XYZ conversion module using RDKit.

This module provides functionality to convert SMILES strings to XYZ format
with explicit hydrogen atoms.
"""

from typing import Final

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles


class ConversionError(Exception):
    """Exception raised when SMILES to XYZ conversion fails."""


def smiles_to_xyz(smiles: str, *, add_hydrogens: bool = True, optimize: bool = True) -> str:
    """Convert a SMILES string to XYZ format.

    Args:
        smiles: The SMILES string to convert.
        add_hydrogens: Whether to add explicit hydrogen atoms. Defaults to True.
        optimize: Whether to optimize the 3D geometry using UFF force field. Defaults to True.

    Returns:
        A string containing the molecule in XYZ format.

    Raises:
        ConversionError: If the SMILES string is invalid or conversion fails.

    Example:
        >>> xyz = smiles_to_xyz("CCO")  # Ethanol
        >>> print(xyz.split("\\n")[0])  # Number of atoms
        9
    """
    # Parse SMILES string
    mol: Chem.Mol | None = Chem.MolFromSmiles(smiles)
    if mol is None:
        msg: Final[str] = f"Invalid SMILES string: {smiles}"
        raise ConversionError(msg)

    # Add explicit hydrogens if requested
    if add_hydrogens:
        mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:  # type: ignore[attr-defined]
        error_msg: str = f"Failed to generate 3D coordinates for SMILES: {smiles}"
        raise ConversionError(error_msg)

    # Optimize geometry using UFF force field
    if optimize and AllChem.UFFOptimizeMolecule(mol) != 0:  # type: ignore[attr-defined]
        error_msg = f"Geometry optimization failed for SMILES: {smiles}"
        raise ConversionError(error_msg)

    # Convert to XYZ format
    return mol_to_xyz(mol)


def mol_to_xyz(mol: Chem.Mol) -> str:
    """Convert an RDKit molecule to XYZ format string.

    Args:
        mol: The RDKit molecule object with 3D coordinates.

    Returns:
        A string containing the molecule in XYZ format.

    Raises:
        ConversionError: If the molecule has no conformer or conversion fails.
    """
    if mol.GetNumConformers() == 0:
        msg: Final[str] = "Molecule has no 3D coordinates (conformer)"
        raise ConversionError(msg)

    # Use RDKit's built-in XYZ block writer
    return rdmolfiles.MolToXYZBlock(mol)
