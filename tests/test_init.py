"""Test package initialization."""

from smiles2xyz import ConversionError, __version__, smiles_to_xyz


def test_version() -> None:
    """Test that version is defined."""
    assert __version__ == "0.1.0"


def test_exports() -> None:
    """Test that main functions are exported."""
    assert callable(smiles_to_xyz)
    assert issubclass(ConversionError, Exception)
