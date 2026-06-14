"""Atom group extensions for ProLint.

This module provides extended atom group classes with additional
functionality for biomolecular interaction analysis.
"""

from abc import ABC, abstractmethod
from typing import Iterable, Union, Dict, Optional, Any
from collections import Counter
from functools import cached_property

import numpy as np
import MDAnalysis as mda


class PLAtomGroupBase(ABC):
    """Abstract base class for ProLint atom group operations.

    Defines the interface for atom group manipulation methods
    used throughout ProLint.
    """

    @abstractmethod
    def add(
        self,
        resname: Optional[Union[str, list[str]]] = None,
        atomname: Optional[Union[str, list[str]]] = None,
        resnum: Optional[Union[int, list[int]]] = None,
        atomids: Optional[Union[int, list[int]]] = None,
    ) -> "ExtendedAtomGroup":
        """Add atoms to the atom group.

        Parameters
        ----------
        resname : str or list of str, optional
            Residue name(s) to add.
        atomname : str or list of str, optional
            Atom name(s) to add.
        resnum : int or list of int, optional
            Residue number(s) to add.
        atomids : int or list of int, optional
            Atom ID(s) to add.

        Returns
        -------
        ExtendedAtomGroup
            New atom group with added atoms.
        """

    @abstractmethod
    def remove(
        self,
        resname: Optional[Union[str, list[str]]] = None,
        atomname: Optional[Union[str, list[str]]] = None,
        resnum: Optional[Union[int, list[int]]] = None,
        atomids: Optional[Union[int, list[int]]] = None,
    ) -> "ExtendedAtomGroup":
        """Remove atoms from the atom group.

        Parameters
        ----------
        resname : str or list of str, optional
            Residue name(s) to remove.
        atomname : str or list of str, optional
            Atom name(s) to remove.
        resnum : int or list of int, optional
            Residue number(s) to remove.
        atomids : int or list of int, optional
            Atom ID(s) to remove.

        Returns
        -------
        ExtendedAtomGroup
            New atom group with specified atoms removed.
        """

    @abstractmethod
    def get_resnames(
        self, resids: Iterable[int], out: Union[type[list], type[dict]] = list
    ) -> Union[list[str], Dict[int, str]]:
        """Get residue names for given residue IDs.

        Parameters
        ----------
        resids : Iterable of int
            Residue IDs to look up.
        out : type, default=list
            Output format: ``list`` or ``dict``.

        Returns
        -------
        list of str or dict
            Residue names as list or as {resid: resname} mapping.
        """

    @abstractmethod
    def get_resids(self, resname: str) -> np.ndarray:
        """Get residue IDs for a given residue name.

        Parameters
        ----------
        resname : str
            Residue name to look up.

        Returns
        -------
        ndarray
            Array of residue IDs matching the residue name.
        """

    @abstractmethod
    def get_all_resids(
        self, resnames: Iterable[str], out: Union[type[list], type[dict]] = list
    ) -> Union[list[np.ndarray], Dict[str, np.ndarray]]:
        """Get residue IDs for multiple residue names.

        Parameters
        ----------
        resnames : Iterable of str
            Residue names to look up.
        out : type, default=list
            Output format: ``list`` or ``dict``.

        Returns
        -------
        list or dict
            Residue IDs as list of arrays or as {resname: resids} mapping.
        """

    @abstractmethod
    def filter_resids_by_resname(
        self, resids: Iterable[int], resname: str
    ) -> np.ndarray:
        """Filter residue IDs to keep only those matching a residue name.

        Parameters
        ----------
        resids : Iterable of int
            Residue IDs to filter.
        resname : str
            Residue name to match.

        Returns
        -------
        ndarray
            Subset of input resids that match the residue name.
        """

    @property
    @abstractmethod
    def unique_resnames(self) -> np.ndarray:
        """Unique residue names in the atom group.

        Returns
        -------
        ndarray
            Array of unique residue names.
        """

    @property
    @abstractmethod
    def resname_counts(self) -> Counter:
        """Count of residues for each residue name.

        Returns
        -------
        Counter
            Mapping of residue name to count.
        """


class ExtendedAtomGroup(mda.AtomGroup, PLAtomGroupBase):
    """Extended MDAnalysis AtomGroup with additional ProLint functionality.

    Provides enhanced methods for manipulating and querying atom selections,
    with cached properties for efficient repeated access.

    Parameters
    ----------
    *args : tuple
        Arguments passed to MDAnalysis AtomGroup.
    **kwargs : dict
        Keyword arguments passed to MDAnalysis AtomGroup.

    Examples
    --------
    >>> from prolint import Universe
    >>> u = Universe("topology.gro", "trajectory.xtc")
    >>> u.database.unique_resnames
    array(['POPC', 'POPE', 'CHOL'], dtype='<U4')
    >>> u.database.resname_counts
    Counter({'POPC': 128, 'POPE': 64, 'CHOL': 32})

    Add atoms:

    >>> extended = u.database.add(resname="DPPC")

    Remove atoms:

    >>> filtered = u.database.remove(resname="CHOL")

    See Also
    --------
    Universe.query : Query atom group (typically protein)
    Universe.database : Database atom group (typically lipids)
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        # Cached properties will be computed on first access
        self._stored_resnames = self.residues.resnames
        self._stored_resids = self.residues.resids

    @cached_property
    def _resname_resid_labels(self) -> Dict[int, str]:
        resnames = self.residues.resnames
        resids = self.residues.resids

        return dict(zip(resids, resnames))

    def _build_selection_string(
        self,
        resname: Optional[Union[str, list[str]]] = None,
        atomname: Optional[Union[str, list[str]]] = None,
        resnum: Optional[Union[int, list[int]]] = None,
        atomids: Optional[Union[int, list[int]]] = None,
    ) -> str:
        """Build MDAnalysis selection string from filter criteria.

        Parameters
        ----------
        resname : str or list of str, optional
            Residue name(s).
        atomname : str or list of str, optional
            Atom name(s).
        resnum : int or list of int, optional
            Residue number(s).
        atomids : int or list of int, optional
            Atom ID(s).

        Returns
        -------
        str
            MDAnalysis selection string.

        Raises
        ------
        ValueError
            If no selection criteria are provided.
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
            if isinstance(resnum, int):
                resnum = [resnum]
            resnum_str = list(map(str, resnum))
            selections.append("resid " + " or resid ".join(resnum_str))

        if atomids is not None:
            if isinstance(atomids, int):
                atomids = [atomids]
            atomids_str = list(map(str, atomids))
            selections.append("bynum " + " or bynum ".join(atomids_str))

        if not selections:
            raise ValueError("At least one selection criterion must be provided")

        return " or ".join(selections)

    def add(
        self,
        resname: Optional[Union[str, list[str]]] = None,
        atomname: Optional[Union[str, list[str]]] = None,
        resnum: Optional[Union[int, list[int]]] = None,
        atomids: Optional[Union[int, list[int]]] = None,
    ) -> "ExtendedAtomGroup":
        """Add atoms to the atom group.

        See :meth:`PLAtomGroupBase.add` for parameter documentation.
        """
        selection_string = self._build_selection_string(
            resname, atomname, resnum, atomids
        )
        new_group = self.universe.atoms.select_atoms(selection_string)
        new_group = self | new_group

        return self.__class__(new_group)

    def remove(
        self,
        resname: Optional[Union[str, list[str]]] = None,
        atomname: Optional[Union[str, list[str]]] = None,
        resnum: Optional[Union[int, list[int]]] = None,
        atomids: Optional[Union[int, list[int]]] = None,
    ) -> "ExtendedAtomGroup":
        """Remove atoms from the atom group.

        See :meth:`PLAtomGroupBase.remove` for parameter documentation.
        """
        selection_string = self._build_selection_string(
            resname, atomname, resnum, atomids
        )
        atoms_to_remove = self.select_atoms(selection_string)
        new_group = self - atoms_to_remove

        return self.__class__(new_group)

    def get_resnames(
        self, resids: Iterable[int], out: Union[type[list], type[dict]] = list
    ) -> Union[list[str], Dict[int, str]]:
        """Get residue names for given residue IDs.

        See :meth:`PLAtomGroupBase.get_resnames` for parameter documentation.
        """
        if out is list:
            return [self._resname_resid_labels[resid] for resid in resids]
        elif out is dict:
            return {resid: self._resname_resid_labels[resid] for resid in resids}
        else:
            raise ValueError("out must be either list or dict")

    def get_resids(self, resname: str) -> np.ndarray:
        """Get residue IDs for a given residue name.

        See :meth:`PLAtomGroupBase.get_resids` for parameter documentation.
        """
        return self.residues.resids[self.residues.resnames == resname]

    def get_all_resids(
        self, resnames: Iterable[str], out: Union[type[list], type[dict]] = list
    ) -> Union[list[np.ndarray], Dict[str, np.ndarray]]:
        """Get residue IDs for multiple residue names.

        See :meth:`PLAtomGroupBase.get_all_resids` for parameter documentation.
        """
        if out is list:
            return [self.get_resids(resname) for resname in resnames]
        elif out is dict:
            return {resname: self.get_resids(resname) for resname in resnames}
        else:
            raise ValueError("out must be either list or dict")

    def filter_resids_by_resname(
        self, resids: Iterable[int], resname: str
    ) -> np.ndarray:
        """Filter residue IDs to keep only those matching a residue name.

        See :meth:`PLAtomGroupBase.filter_resids_by_resname` for parameter documentation.
        """
        resids = np.asarray(resids)
        all_resnames = self._stored_resnames
        all_resids = self._stored_resids
        indices = np.searchsorted(all_resids, resids)
        return resids[np.where(all_resnames[indices] == resname)[0]]

    @property
    def unique_resnames(self) -> np.ndarray:
        """Unique residue names in the atom group.

        Returns
        -------
        ndarray
            Array of unique residue names.
        """
        return np.unique(self.residues.resnames)  # type: ignore[return-value]

    @property
    def resname_counts(self) -> Counter:
        """Count of residues for each residue name.

        Returns
        -------
        Counter
            Mapping of residue name to count.
        """
        return Counter(self.residues.resnames)

    def __str__(self) -> str:
        return f"<ProLint Wrapper for {super().__str__()}>"

    def __repr__(self) -> str:
        return f"<ProLint Wrapper for {super().__repr__()}>"
