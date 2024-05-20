r""":mod:`prolint2.core.contact_provider`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from collections import defaultdict
from typing import Callable, Literal

import numpy as np
import pandas as pd

from prolint2.computers.contacts import ContactComputerBase, SerialContacts
from prolint2.core.typing import (
    NestedFloatDict,
    NestedIterFloatDict,
    NestedIterIntDict,
    LipidId,
)

from prolint2.metrics.base import BaseContactStore
from prolint2.metrics.exact_contacts import ExactContacts
from prolint2.metrics.aprox_contacts import AproxContacts

from prolint2.config.units import DEFAULT_SIM_PARAMS


class ComputedContacts:
    """A class to compute contacts between residues and lipids.

    Parameters
    ----------
    contact_strategy_instance : BaseContactStore
        An instance of a contact strategy class.
    provider : ContactsProvider
        The contact provider that will be used to compute contacts.

    """

    def __init__(
        self, contact_strategy_instance: BaseContactStore, provider: "ContactsProvider"
    ):
        self._contact_strategy = contact_strategy_instance
        self.provider = provider

    def compute_metric(self, metric: str, target_lipid_name=None) -> NestedFloatDict:
        """Compute a pre-defined metric for all lipids or a specific lipid.

        Parameters
        ----------
        metric : str
            The metric to compute. Must be one of 'max', 'sum', 'mean'.
        target_lipid_name : str, optional
            The name of the lipid to compute the metric for. If None, the metric will be computed for all lipids.

        Returns
        -------
        Dict[str, Dict[str, Dict[int, float]]]
            A dictionary of computed metrics for all lipids.

        Examples
        --------
        >>> c.compute('max')
        >>> c.compute('sum', 'DOPC')
        >>> c.compute('median') # raises ValueError. Use `apply_function` instead.
        """

        return self._contact_strategy.compute(
            metric, target_lipid_name=target_lipid_name
        )

    def apply_function(self, func: Callable, target_lipid_name=None) -> NestedFloatDict:
        """Apply the given function to the contacts for the given lipid name."""
        return self._contact_strategy.apply_function(
            func, target_lipid_name=target_lipid_name
        )

    @property
    def contacts(self) -> NestedIterFloatDict:
        """The computed contacts."""
        return self._contact_strategy.contacts

    @property
    def pooled_contacts(self) -> NestedIterFloatDict:
        """The computed contacts."""
        return self._contact_strategy.pooled_results()

    @property
    def contact_frames(self) -> NestedIterIntDict:
        """The computed contacts."""
        return self._contact_strategy.contact_frames
    
    @property
    def occupancies(self) -> NestedIterIntDict:
        """The computed contacts."""
        return self._contact_strategy.occupancies

    def create_dataframe(self, n_frames: int) -> pd.DataFrame:
        """Create a pandas DataFrame from the computed contacts.

        Parameters
        ----------
        n_frames : int
            The number of frames in the trajectory.

        Returns
        -------
        pd.DataFrame
            A pandas DataFrame with the computed contacts.
        """
        keys = []
        contact_arrays = []

        for residue_id, lipid_name_dict in self.contact_frames.items():
            for lipid_id, frame_indices in lipid_name_dict.items():
                contact_array = np.zeros(n_frames, dtype=np.int8)
                contact_array[frame_indices] = 1

                keys.append((residue_id, lipid_id))
                contact_arrays.append(contact_array)

        df = pd.DataFrame(
            contact_arrays,
            index=pd.MultiIndex.from_tuples(keys, names=["ResidueID", "LipidId"]),
        )
        df = df.sort_index(level=["ResidueID", "LipidId"], ascending=[True, True])

        return df

    def get_lipids_by_residue_id(self, residue_id: int) -> list:
        """Get all LipidIds that interact with the given ResidueID."""
        return sorted(list(self.contact_frames[residue_id].keys()))

    def get_residues_by_lipid_id(self, lipid_id: int) -> list:
        """Get all ResidueIDs that interact with the given LipidId."""
        residues = [
            residue_id
            for residue_id, lipid_name_dict in self.contact_frames.items()
            if lipid_id in lipid_name_dict.keys()
        ]
        return residues

    def get_contact_data(
        self, residue_id: int, lipid_id: int, output: str = "contacts"
    ) -> list:
        """Get the contact data for a given residue and lipid.

        Parameters
        ----------
        residue_id : int
            The residue id.
        lipid_id : int
            The lipid id.
        output : str, optional
            The output format. Must be one of 'contacts' or 'indices'.

        Returns
        -------
        list
            A list of contacts or frame indices.
        """

        frame_indices = self.contact_frames[residue_id][lipid_id]

        if output == "indices":
            return frame_indices
        else:
            n_frames = (
                max(
                    [
                        max(frame_indices_list)
                        for frame_indices_list in self.contact_frames[
                            residue_id
                        ].values()
                    ]
                )
                + 1
            )
            contact_array = [1 if i in frame_indices else 0 for i in range(n_frames)]
            return contact_array

    def intersection(self, other: "ComputedContacts") -> "ComputedContacts":
        """Compute the intersection of two contact providers. Note that ProLint contacts use a radial cutoff.
        This means that the intersection between two contact providers (c1 and c2) will be equal to the contact provider
        with the smallest cutoff. ProLint, however, defines the intersection between two contact providers (c1 and c2) to
        be equal to the lipid ids of the contact provider with the smallest cutoff, and the frame indices of the contact
        provider with the largest cutoff. This way the intersection between two contact providers is meaningful and
        computationaly allows for chaining of contact providers (See example below).

        Parameters
        ----------
        other : ComputedContacts
            The other contact provider to compute the intersection with.

        Returns
        -------
        ContactsProvider
            A new contact provider with the intersection of the contacts of both contact providers.

        Examples
        --------
        >>> ts = Universe('coordinates.gro', 'trajectory.xtc')
        >>> c1 = ts.compute_contacts(cutoff=7)
        >>> c2 = ts.compute_contacts(cutoff=8)
        >>> c3 = c1 + c2
        >>> c1 + c2 == c2 + c1 # True
        """
        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, lipid_ids in self.contact_frames.items():
            for lipid_id in lipid_ids:
                if LipidId(lipid_id) in other.contact_frames[residue_id]:
                    result_data[residue_id][lipid_id] = other.contact_frames[
                        residue_id
                    ][lipid_id]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(
            self.provider.query.universe, result_data
        )
        contact_instances.norm_factor = self.provider.params.get("norm_factor", 1)
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def difference(self, other: "ComputedContacts") -> "ComputedContacts":
        """Compute the difference of two contact providers. Given two contact providers (c1 and c2), the difference
        between them (c2 -c1) is defined as the contacts of c2 that are not present in c1.

        Parameters
        ----------
        other : ComputedContacts
            The other contact provider to compute the difference with.

        Returns
        -------
        ContactsProvider
            A new contact provider with the difference of the contacts of both contact providers.

        Examples
        --------
        >>> ts = Universe('coordinates.gro', 'trajectory.xtc')
        >>> c1 = ts.compute_contacts(cutoff=7)
        >>> c2 = ts.compute_contacts(cutoff=8)
        >>> c3 = c2 - c1
        >>> c1 - c2 == c2 - c1 # False, c1 - c2 will be an empty contact provider if c1 is a subset of c2
        """

        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, lipid_ids in self.contact_frames.items():
            for lipid_id in lipid_ids:
                if LipidId(lipid_id) not in other.contact_frames[residue_id]:
                    result_data[residue_id][lipid_id] = self.contact_frames[residue_id][
                        lipid_id
                    ]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(
            self.provider.query.universe, result_data
        )
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def __add__(self, other: "ComputedContacts") -> "ComputedContacts":
        return self.intersection(other)

    def __sub__(self, other: "ComputedContacts") -> "ComputedContacts":
        return self.difference(other)


class ContactsProvider:
    """
    Class that provides the contacts computation functionality.
    """

    def __init__(
        self,
        query,
        database,
        params=None,
        compute_strategy: Literal["default"] = "default",
        contact_strategy: Literal["exact", "aprox"] = "exact",
    ):
        self.query = query
        self.database = database

        self._contact_computers = {"default": SerialContacts}
        self._contact_counter = {"exact": ExactContacts, "aprox": AproxContacts}
        self._compute_strategy = compute_strategy
        self._contact_strategy = self._contact_counter[contact_strategy]

        self.params = params if params is not None else DEFAULT_SIM_PARAMS

    def compute(
        self, strategy_or_computer=None, start=None, stop=None, step=1, **kwargs
    ):
        """
        Compute contacts between the query and the database.

        Parameters
        ----------
        strategy_or_computer : str or ContactComputerBase, optional
            The strategy to compute contacts. If None, the default strategy is used.
        **kwargs
            Additional arguments to pass to the contact computer.

        Returns
        -------
        ComputedContacts
            The computed contacts.
        """
        if strategy_or_computer is None:
            strategy_or_computer = self._compute_strategy

        # Strategy to compute contacts (e.g. serial, parallel, etc.)
        if isinstance(strategy_or_computer, ContactComputerBase):
            contact_computer = strategy_or_computer
        else:
            contact_computer_class = self._contact_computers.get(
                strategy_or_computer, None
            )
            if contact_computer_class is None:
                strats = ", ".join(self._contact_computers.keys())
                raise ValueError(
                    f"Unknown strategy or computer: {strategy_or_computer}. Available strategies are: {strats}."
                )
            contact_computer = contact_computer_class(
                self.query.universe, self.query, self.database, **kwargs
            )
        contact_computer.run(verbose=True, start=start, stop=stop, step=step)

        # Strategy to count and store contacts (e.g. exact, aprox, etc.)
        contact_strategy_instance = self._contact_strategy(
            self.query.universe,
            contact_computer.contact_frames,
            self.params.get("norm_factor"),
        )
        contact_strategy_instance.run()

        return ComputedContacts(contact_strategy_instance, self)

    def load_from_file(self, file, **kwargs):
        """
        Load contacts from a file.

        Parameters
        ----------
        file : str or pathlib.Path
            The path to the file to load the contacts from.
        **kwargs
            Additional arguments to pass to the contact loader.

        Returns
        -------
        ComputedContacts
            The computed contacts.
        """
        # get contact frames from file
        df = pd.read_csv(file, index_col=[0, 1])
        contact_frames = defaultdict(lambda: defaultdict(list))
        for residue_id, lipid_id in df.index:
            contact_frames[residue_id][lipid_id] = np.nonzero(
                df.loc[(residue_id, lipid_id)].to_numpy()
            )[0].tolist()

        # Count and store contacts
        contact_strategy_instance = self._contact_strategy(
            self.query.universe, contact_frames, self.params.get("norm_factor")
        )
        contact_strategy_instance.run()

        return ComputedContacts(contact_strategy_instance, self)
