"""Tests for CLI functionality."""

import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

from smiles2xyz.cli import create_parser, main


class TestCLI:
    """Test suite for command-line interface."""

    def test_create_parser(self) -> None:
        """Test that parser is created correctly."""
        parser = create_parser()
        assert parser.prog == "smiles2xyz"

    def test_parser_required_argument(self) -> None:
        """Test that SMILES is a required argument."""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_basic_usage(self) -> None:
        """Test basic parser usage."""
        parser = create_parser()
        args = parser.parse_args(["CCO"])
        assert args.smiles == "CCO"
        assert args.hydrogens is True
        assert args.optimize is True
        assert args.output is None

    def test_parser_with_output(self) -> None:
        """Test parser with output file."""
        parser = create_parser()
        args = parser.parse_args(["CCO", "-o", "output.xyz"])
        assert args.output == "output.xyz"

    def test_parser_no_hydrogens(self) -> None:
        """Test parser with --no-hydrogens flag."""
        parser = create_parser()
        args = parser.parse_args(["CCO", "--no-hydrogens"])
        assert args.hydrogens is False

    def test_parser_no_optimize(self) -> None:
        """Test parser with --no-optimize flag."""
        parser = create_parser()
        args = parser.parse_args(["CCO", "--no-optimize"])
        assert args.optimize is False

    def test_main_stdout(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test main function with stdout output."""
        test_args = ["smiles2xyz", "C"]
        with patch.object(sys, "argv", test_args):
            main()

        captured = capsys.readouterr()
        # Should output XYZ format
        assert "5" in captured.out  # 5 atoms for CH4
        assert "C " in captured.out

    def test_main_with_file_output(self, tmp_path: Path) -> None:
        """Test main function with file output."""
        output_file = tmp_path / "test.xyz"
        test_args = ["smiles2xyz", "C", "-o", str(output_file)]

        with patch.object(sys, "argv", test_args):
            main()

        assert output_file.exists()
        content = output_file.read_text()
        assert "5" in content  # 5 atoms for CH4

    def test_main_invalid_smiles(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test main function with invalid SMILES."""
        test_args = ["smiles2xyz", "INVALID!!!"]

        with patch.object(sys, "argv", test_args):
            with pytest.raises(SystemExit) as exc_info:
                main()

        assert exc_info.value.code == 1
        captured = capsys.readouterr()
        assert "Error:" in captured.err
