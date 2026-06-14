"""Unit tests for prolint.computers.base module."""

from abc import ABC


class TestContactComputerBase:
    """Tests for ContactComputerBase class."""

    def test_inherits_from_analysis_base_and_abc(self):
        """ContactComputerBase inherits from AnalysisBase and ABC."""
        from prolint.computers.base import ContactComputerBase
        from MDAnalysis.analysis.base import AnalysisBase

        assert issubclass(ContactComputerBase, AnalysisBase)
        assert issubclass(ContactComputerBase, ABC)
