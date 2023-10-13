r""":mod:`prolint2.core.base`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import numpy as np
from MDAnalysis.core.topologyattrs import ResidueStringAttr


class MacrosClass(ResidueStringAttr):
    """
    A class for managing macro attributes associated with residues in a universe.

    This class inherits from ResidueStringAttr and is designed for managing macro attributes.

    :param universe: The universe to which this class is associated.
    :type universe: Universe
    """

    attrname = "macros"
    singular = "macro"

    def __init__(self, universe):
        """
        Initialize MacrosClass with the given universe.

        :param universe: The universe to which this class is associated.
        :type universe: Universe
        """
        n_atoms, n_residues, n_segments = (
            universe.atoms.n_atoms,
            universe.residues.n_residues,
            universe.segments.n_segments,
        )
        values = self._gen_initial_values(n_atoms, n_residues, n_segments)
        super().__init__(values)

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments):
        """
        Generate initial macro attribute values.

        :param n_atoms: The number of atoms in the universe.
        :type n_atoms: int
        :param n_residues: The number of residues in the universe.
        :type n_residues: int
        :param n_segments: The number of segments in the universe.
        :type n_segments: int
        :return: An array of initial macro attribute values.
        :rtype: numpy.ndarray
        """
        return np.array(["other"] * n_residues, dtype=object)

    @staticmethod
    def set_macros_values(query, n_proteins=None):
        """
        Set macro attribute values for residues based on the query.

        :param query: The query used to set macro attribute values.
        :type query: Universe
        :param n_proteins: The number of proteins (optional).
        :type n_proteins: int
        """
        protein_segments = query.segments

        if len(protein_segments) == 1 and protein_segments.n_atoms == query.n_atoms:
            for segment_idx, segment in enumerate(protein_segments):
                segment.residues.macros = "Protein" + str(segment_idx)
        else:
            resseq = query.residues.resids
            res0, first_last_index, first_index = resseq[0], [], 0

            for last_index, res in enumerate(resseq):
                if res < res0:
                    first_last_index.append((first_index, last_index - 1))
                    first_index = last_index
                res0 = res
            first_last_index.append((first_index, last_index))

            for idx, (first_index, last_index) in enumerate(first_last_index):
                selected_residues = query.residues[first_index : last_index + 1]
                selected_residues.macros = "Protein" + str(idx)
