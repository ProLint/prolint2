r""":mod:`prolint2.core.groups`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from abc import ABC, abstractmethod
from typing import Iterable, Union, Dict
from collections import Counter

import numpy as np
import MDAnalysis as mda


class PLAtomGroupBase(ABC):
    """
    An abstract base class for AtomGroup objects.

    This abstract base class defines a set of methods and properties for working with AtomGroup objects.

    """

    @abstractmethod
    def add(self, resname=None, atomname=None, resnum=None, atomids=None):
        """
        Add atoms to the query or database.

        :param resname: The residue name to filter by.
        :param atomname: The atom name to filter by.
        :param resnum: The residue number to filter by.
        :param atomids: A list of atom IDs to add.
        :type resname: str, optional
        :type atomname: str, optional
        :type resnum: int, optional
        :type atomids: List[int], optional

        """

    @abstractmethod
    def remove(self, resname=None, atomname=None, resnum=None, atomids=None):
        """
        Remove atoms from the query or database.

        :param resname: The residue name to filter by.
        :param atomname: The atom name to filter by.
        :param resnum: The residue number to filter by.
        :param atomids: A list of atom IDs to remove.
        :type resname: str, optional
        :type atomname: str, optional
        :type resnum: int, optional
        :type atomids: List[int], optional

        """

    @abstractmethod
    def get_resname(self, resid: int):
        """
        Get the residue name of a residue in the AtomGroup.

        :param resid: The residue ID.
        :type resid: int
        :return: The residue name of the specified residue.
        :rtype: str

        """

    @abstractmethod
    def get_resnames(self, resids: Iterable[int]):
        """
        Get the residue names of a list of residues in the AtomGroup.

        :param resids: An iterable of residue IDs.
        :type resids: Iterable[int]
        :return: A list of residue names corresponding to the input residue IDs.
        :rtype: List[str]

        """

    @abstractmethod
    def get_resid(self, resname: str):
        """
        Get the residue ID of a residue in the AtomGroup.

        :param resname: The residue name.
        :type resname: str
        :return: The residue ID of the specified residue name.
        :rtype: int

        """

    @abstractmethod
    def get_resids(self, resnames: Iterable[str]):
        """
        Get the residue IDs of a list of residues in the AtomGroup.

        :param resnames: An iterable of residue names.
        :type resnames: Iterable[str]
        :return: A list of residue IDs corresponding to the input residue names.
        :rtype: List[int]

        """

    @abstractmethod
    def filter_resids_by_resname(self, resids: np.ndarray, resname: str):
        """
        Filter the residue IDs by residue name.

        :param resids: An array of residue IDs to filter.
        :param resname: The residue name to filter by.
        :type resids: np.ndarray
        :type resname: str
        :return: An array of residue IDs that match the specified residue name.
        :rtype: np.ndarray

        """

    @property
    @abstractmethod
    def unique_resnames(self):
        """
        Get the unique residue names in the AtomGroup.

        :return: A set of unique residue names present in the AtomGroup.
        :rtype: set

        """

    @property
    @abstractmethod
    def resname_counts(self):
        """
        Get the number of residues of each residue name in the AtomGroup.

        :return: A dictionary containing the counts of each unique residue name.
        :rtype: dict

        """

        """ Get the number of residues of each residue name in the AtomGroup."""


class ExtendedAtomGroup(mda.AtomGroup, PLAtomGroupBase):
    """
    An extended version of the MDAnalysis AtomGroup class.

    This class extends the functionality of the MDAnalysis AtomGroup to provide additional
    features for working with atoms and residues.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the AtomGroup.

        Parameters:
        *args: Positional arguments for the superclass constructor.
        **kwargs: Keyword arguments for the superclass constructor.
        """
        super().__init__(*args, **kwargs)
        self._resname_resid_labels = self._build_resname_resid_labels()
        self._stored_resnames = self.residues.resnames
        self._stored_resids = self.residues.resids

    def _build_resname_resid_labels(self):
        """
        Build a dictionary of residue names and residue IDs.

        Returns:
        dict: A dictionary mapping residue IDs to residue names.
        """
        resnames = self.residues.resnames
        resids = self.residues.resids

        return dict(zip(resids, resnames))

    def _build_stored_resnames(self):
        """
        Build a dictionary of residue names and residue IDs.

        Returns:
        list: A list of residue names.
        """
        resnames = self.residues.resnames
        return resnames

    def _build_selection_string(
        self, resname=None, atomname=None, resnum=None, atomids=None
    ):
        """
        Build a selection string for atom filtering.

        Parameters:
        resname (str or list): Residue name(s) for filtering.
        atomname (str or list): Atom name(s) for filtering.
        resnum (int or list): Residue number(s) for filtering.
        atomids (int or list): Atom IDs for filtering.

        Returns:
        str: A selection string for atom filtering.
        Raises:
        ValueError: If no selection criteria are provided.
        """
        selections = []

        if resname is not None:
            if isinstance(resname, str):
                resname = [resname]
            selections.append("resname " + " or resname ".join(resname))

        if atomname is not None:
            if isinstance(atomname, str):
                atomname = [atomname]
            selections.append("name " + " or name ".join(atomname))

        if resnum is not None:
            resnum = map(str, resnum)
            selections.append("resid " + " or resid ".join(resnum))

        if atomids is not None:
            atomids = map(str, atomids)
            selections.append("bynum " + " or bynum ".join(atomids))

        if not selections:
            raise ValueError("At least one selection criterion must be provided")

        return " or ".join(selections)

    def add(self, resname=None, atomname=None, resnum=None, atomids=None):
        """
        Add atoms to the query or database.

        Parameters:
        resname (str or list): Residue name(s) for filtering.
        atomname (str or list): Atom name(s) for filtering.
        resnum (int or list): Residue number(s) for filtering.
        atomids (int or list): Atom IDs for filtering.

        Returns:
        ExtendedAtomGroup: A new ExtendedAtomGroup containing the selected atoms.
        """
        selection_string = self._build_selection_string(
            resname, atomname, resnum, atomids
        )
        new_group = self.universe.atoms.select_atoms(selection_string)
        new_group = self | new_group

        return self.__class__(new_group)

    def remove(self, resname=None, atomname=None, resnum=None, atomids=None):
        """
        Remove atoms from the query or database.

        Parameters:
        resname (str or list): Residue name(s) for filtering.
        atomname (str or list): Atom name(s) for filtering.
        resnum (int or list): Residue number(s) for filtering.
        atomids (int or list): Atom IDs for filtering.

        Returns:
        ExtendedAtomGroup: A new ExtendedAtomGroup with the selected atoms removed.
        """
        selection_string = self._build_selection_string(
            resname, atomname, resnum, atomids
        )
        atoms_to_remove = self.select_atoms(selection_string)
        new_group = self - atoms_to_remove

        return self.__class__(new_group)

    def get_resname(self, resid: int):
        """
        Get the residue name of a residue in the AtomGroup.

        Parameters:
        resid (int): Residue ID.

        Returns:
        str: Residue name associated with the given residue ID.
        """
        return self._resname_resid_labels[resid]

    def get_resnames(
        self, resids: Iterable[int], out: Union[list, Dict[int, str]] = list
    ):
        """
        Get the residue names of a list of residues in the AtomGroup.

        Parameters:
        resids (iterable): List of residue IDs.
        out (list or dict): Output format.

        Returns:
        list or dict: Residue names corresponding to the given residue IDs.
        Raises:
        ValueError: If the 'out' parameter is not 'list' or 'dict'.
        """
        if out is list:
            return [self._resname_resid_labels[resid] for resid in resids]
        elif out is dict:
            return {resid: self._resname_resid_labels[resid] for resid in resids}
        else:
            raise ValueError("out must be either list or dict")

    def get_resid(self, resname: str):
        """
        Get the residue ID of a residue in the AtomGroup.

        Parameters:
        resname (str): Residue name.

        Returns:
        int: Residue ID associated with the given residue name.
        """
        return self.residues.resids[self.residues.resnames == resname][0]

    def get_resids(
        self, resnames: Iterable[str], out: Union[list, Dict[str, int]] = list
    ):
        """
        Get the residue IDs of a list of residues in the AtomGroup.

        Parameters:
        resnames (iterable): List of residue names.
        out (list or dict): Output format.

        Returns:
        list or dict: Residue IDs corresponding to the given residue names.
        Raises:
        ValueError: If the 'out' parameter is not 'list' or 'dict'.
        """
        if out is list:
            return [self.get_resid(resname) for resname in resnames]
        elif out is dict:
            return {resname: self.get_resid(resname) for resname in resnames}
        else:
            raise ValueError("out must be either list or dict")

    def filter_resids_by_resname(self, resids: Iterable[int], resname: str):
        """
        Filter the residue IDs by residue name.

        Parameters:
        resids (iterable): List of residue IDs to filter.
        resname (str): Residue name for filtering.

        Returns:
        np.ndarray: Filtered residue IDs.
        """
        resids = np.asarray(resids)
        all_resnames = self._stored_resnames
        all_resids = self._stored_resids
        # print ('shapes', all_resnames.shape, all_resids.shape, resids.shape)
        indices = np.searchsorted(all_resids, resids)
        return resids[np.where(all_resnames[indices] == resname)[0]]

    @staticmethod
    def static_filter_resids_by_resname(
        resids: np.ndarray,
        resnames: np.ndarray,
        resids_subset: np.ndarray,
        resname: str,
    ):
        """
        Filter the residue IDs by residue name using a static method.

        Parameters:
        resids (np.ndarray): All residue IDs in the AtomGroup.
        resnames (np.ndarray): All residue names in the AtomGroup.
        resids_subset (np.ndarray): Subset of residue IDs to filter.
        resname (str): Residue name for filtering.

        Returns:
        np.ndarray: Filtered residue IDs from the subset.
        """
        indices = np.searchsorted(resids, resids_subset)
        return resids_subset[np.where(resnames[indices] == resname)[0]]

    @property
    def unique_resnames(self):
        """
        Get the unique residue names in the AtomGroup.

        Returns:
        np.ndarray: Array of unique residue names.
        """
        return np.unique(self.residues.resnames)

    @property
    def resname_counts(self):
        """
        Get the number of residues of each residue name in the AtomGroup.

        Returns:
        collections.Counter: Counter object with residue name counts.
        """
        return Counter(self.residues.resnames)

    def __str__(self) -> str:
        """
        Return a string representation of the ExtendedAtomGroup.

        Returns:
        str: String representation of the object.
        """
        return f"<ProLint Wrapper for {super().__str__()}>"

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the ExtendedAtomGroup.

        Returns:
        str: Detailed string representation of the object.
        """
        return f"<ProLint Wrapper for {super().__repr__()}>"
