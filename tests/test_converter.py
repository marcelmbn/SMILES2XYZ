"""Tests for the SMILES to XYZ converter."""

import pytest
from rdkit import Chem

from smiles2xyz import ConversionError, smiles_to_xyz
from smiles2xyz.converter import mol_to_xyz


class TestSMILESToXYZ:
    """Test suite for SMILES to XYZ conversion."""

    def test_simple_molecule_ethanol(self) -> None:
        """Test conversion of ethanol (CCO)."""
        xyz = smiles_to_xyz("CCO")
        lines = xyz.split("\n")

        # Check number of atoms (2 C + 1 O + 6 H = 9 atoms)
        assert lines[0] == "9"
        # Second line is comment (empty)
        assert lines[1] == ""
        # Should have 9 atom coordinate lines
        assert len(lines) == 11  # num_atoms + comment + 9 atom lines

    def test_simple_molecule_water(self) -> None:
        """Test conversion of water (O)."""
        xyz = smiles_to_xyz("O")
        lines = xyz.split("\n")

        # Water with hydrogens: 1 O + 2 H = 3 atoms
        assert lines[0] == "3"
        assert lines[1] == ""
        assert len(lines) == 5  # num_atoms + comment + 3 atom lines

    def test_simple_molecule_methane(self) -> None:
        """Test conversion of methane (C)."""
        xyz = smiles_to_xyz("C")
        lines = xyz.split("\n")

        # Methane: 1 C + 4 H = 5 atoms
        assert lines[0] == "5"
        assert "C " in lines[2]  # First atom should be carbon

    def test_benzene(self) -> None:
        """Test conversion of benzene (c1ccccc1)."""
        xyz = smiles_to_xyz("c1ccccc1")
        lines = xyz.split("\n")

        # Benzene: 6 C + 6 H = 12 atoms
        assert lines[0] == "12"

    def test_without_hydrogens(self) -> None:
        """Test conversion without adding explicit hydrogens."""
        xyz = smiles_to_xyz("CCO", add_hydrogens=False)
        lines = xyz.split("\n")

        # Without hydrogens: 2 C + 1 O = 3 atoms
        assert lines[0] == "3"

    def test_without_optimization(self) -> None:
        """Test conversion without geometry optimization."""
        xyz = smiles_to_xyz("CCO", optimize=False)
        lines = xyz.split("\n")

        # Should still have correct number of atoms
        assert lines[0] == "9"

    def test_invalid_smiles(self) -> None:
        """Test that invalid SMILES raises ConversionError."""
        with pytest.raises(ConversionError, match="Invalid SMILES string"):
            smiles_to_xyz("not_a_valid_smiles!!!")

    def test_xyz_format(self) -> None:
        """Test that output is in valid XYZ format."""
        xyz = smiles_to_xyz("C")
        lines = xyz.split("\n")

        # Check first line is a number
        assert lines[0].isdigit()

        # Check coordinate lines have correct format
        for line in lines[2:]:
            if line.strip():  # Skip empty lines
                parts = line.split()
                assert len(parts) == 4  # symbol + x + y + z
                # First part should be an element symbol
                assert parts[0] in ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"]
                # Rest should be floatable
                float(parts[1])
                float(parts[2])
                float(parts[3])

    def test_coordinates_are_3d(self) -> None:
        """Test that generated coordinates are 3D (not all zero)."""
        xyz = smiles_to_xyz("CCO")
        lines = xyz.split("\n")

        # Extract coordinates
        coords = []
        for line in lines[2:]:
            if line.strip():
                parts = line.split()
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        # Check that not all coordinates are (0, 0, 0)
        non_zero_coords = [c for c in coords if c != (0.0, 0.0, 0.0)]
        assert len(non_zero_coords) > 0


class TestMolToXYZ:
    """Test suite for RDKit Mol to XYZ conversion."""

    def test_mol_without_conformer_raises_error(self) -> None:
        """Test that molecule without 3D coords raises ConversionError."""
        mol = Chem.MolFromSmiles("CCO")
        with pytest.raises(ConversionError, match="no 3D coordinates"):
            mol_to_xyz(mol)

    def test_mol_with_hydrogens(self) -> None:
        """Test conversion of molecule with explicit hydrogens."""
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        from rdkit.Chem import AllChem

        AllChem.EmbedMolecule(mol, randomSeed=42)

        xyz = mol_to_xyz(mol)
        lines = xyz.split("\n")

        assert lines[0] == "5"  # CH4 has 5 atoms
