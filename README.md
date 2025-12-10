# SMILES2XYZ

Python package for CLI conversion of SMILES to XYZ formats with explicit hydrogen atoms.

## Features

- Convert SMILES strings to XYZ coordinate format
- Automatic 3D structure generation with RDKit
- Explicit hydrogen atoms added by default
- UFF force field geometry optimization
- Command-line interface and Python API
- Type-safe implementation with full type hints (Python 3.12+)

## Installation

```bash
pip install -e .
```

## Requirements

- Python >= 3.12
- RDKit >= 2025.09.3

## Usage

### Command Line

```bash
# Basic usage (output to stdout)
smiles2xyz "CCO"

# Save to file
smiles2xyz "CCO" -o ethanol.xyz

# Without explicit hydrogens
smiles2xyz "CCO" --no-hydrogens -o ethanol_no_h.xyz

# Skip geometry optimization
smiles2xyz "CCO" --no-optimize -o ethanol_unopt.xyz
```

### Python API

```python
from smiles2xyz import smiles_to_xyz, ConversionError

# Convert SMILES to XYZ format
try:
    xyz_content = smiles_to_xyz("CCO")  # Ethanol
    print(xyz_content)
except ConversionError as e:
    print(f"Conversion failed: {e}")

# Without hydrogens or optimization
xyz_content = smiles_to_xyz("CCO", add_hydrogens=False, optimize=False)
```

## Development

```bash
# Install with development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Type checking
mypy smiles2xyz

# Linting
ruff check smiles2xyz
```

## License

GPL-3.0 - see LICENSE file for details.
