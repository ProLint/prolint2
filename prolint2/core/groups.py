from abc import ABC, abstractmethod
from typing import Iterable, Union, Dict

import MDAnalysis as mda

class PLAtomGroupBase(ABC):
    """An abstract base class for AtomGroup objects."""

    @abstractmethod
    def add(self, resname=None, atomname=None, resnum=None, atomids=None):
        """ Add atoms to the query or database."""

    @abstractmethod
    def remove(self, resname=None, atomname=None, resnum=None, atomids=None):
        """ Remove atoms from the query or database."""

    @abstractmethod
    def get_resname(self, resid: int):
        """ Get the residue name of a residue in the AtomGroup."""

    @abstractmethod
    def get_resnames(self, resids: Iterable[int]):
        """ Get the residue names of a list of residues in the AtomGroup."""

class ExtendedAtomGroup(mda.AtomGroup, PLAtomGroupBase):
    """An extended version of the MDAnalysis AtomGroup class."""

    def __init__(self, *args, **kwargs):
        """Initialize the AtomGroup."""
        super().__init__(*args, **kwargs)
        self._resname_resid_labels = self._build_resname_resid_labels()

    def _build_resname_resid_labels(self):
        """Build a dictionary of residue names and residue IDs."""
        resnames = self.residues.resnames
        resids = self.residues.resids

        return dict(zip(resids, resnames))

    def _build_selection_string(self, resname=None, atomname=None, resnum=None, atomids=None):
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
        """Add atoms to the query or database."""
        selection_string = self._build_selection_string(resname, atomname, resnum, atomids)
        new_group = self.select_atoms(selection_string)
        new_group = self | new_group

        return self.__class__(new_group)

    def remove(self, resname=None, atomname=None, resnum=None, atomids=None):
        """Remove atoms from the query or database."""
        selection_string = self._build_selection_string(resname, atomname, resnum, atomids)
        atoms_to_remove = self.select_atoms(selection_string)
        new_group = self - atoms_to_remove

        return self.__class__(new_group)

    def get_resname(self, resid: int):
        """Get the residue name of a residue in the AtomGroup."""
        return self._resname_resid_labels[resid]

    def get_resnames(self, resids: Iterable[int], out: Union[list, Dict[int, str]] = list):
        """Get the residue names of a list of residues in the AtomGroup."""
        if out is list:
            return [self._resname_resid_labels[resid] for resid in resids]
        elif out is dict:
            return {resid: self._resname_resid_labels[resid] for resid in resids}
        else:
            raise ValueError("out must be either list or dict")

    def __str__(self):
        """Return a string representation of the AtomGroup."""
        return f"ProLint AtomGroup with {len(self)} atoms"

    def __repr__(self):
        """Return a string representation of the AtomGroup."""
        return f"ProLint AtomGroup with {len(self)} atoms"
