from copy import deepcopy
from collections import defaultdict
from typing import Callable, Literal, Union
from prolint2.computers.contacts import ContactComputerBase, SerialContacts
from prolint2.core.typing import NestedFloatDict, NestedIterFloatDict, NestedIterIntDict, LipidId

from prolint2.metrics.exact_contacts import ExactContacts
from prolint2.metrics.aprox_contacts import AproxContacts

from prolint2.config.units import DEFAULT_SIM_PARAMS


class ComputedContacts:
    def __init__(self, contact_strategy_instance: Union[ExactContacts, AproxContacts], provider: 'ContactsProvider'):
        self._contact_strategy = contact_strategy_instance
        self.provider = provider

    def compute_metric(self, metric: str, target_lipid_name=None) -> NestedFloatDict:
        return self._contact_strategy.compute(metric, target_lipid_name=target_lipid_name)

    def apply_function(self, func: Callable, target_lipid_name=None) -> NestedFloatDict:
        return self._contact_strategy.apply_function(func, target_lipid_name=target_lipid_name)

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

    def intersection(self, other: 'ComputedContacts') -> 'ComputedContacts':
        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, lipid_ids in self.contact_frames.items():
            for lipid_id in lipid_ids:
                if LipidId(lipid_id) in other.contact_frames[residue_id]:
                    result_data[residue_id][lipid_id] = other.contact_frames[residue_id][lipid_id]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(self.provider.query.universe, deepcopy(result_data))
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def difference(self, other: 'ComputedContacts') -> 'ComputedContacts':
        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, lipid_ids in self.contact_frames.items():
            for lipid_id in lipid_ids:
                if LipidId(lipid_id) not in other.contact_frames[residue_id]:
                    result_data[residue_id][lipid_id] = self.contact_frames[residue_id][lipid_id]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(self.provider.query.universe, deepcopy(result_data))
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def __add__(self, other: 'ComputedContacts') -> 'ComputedContacts':
        return self.intersection(other)
    
    def __sub__(self, other: 'ComputedContacts') -> 'ComputedContacts':
        return self.difference(other)


class ContactsProvider:
    """
    Class that provides the contacts computation functionality.
    """
    def __init__(self, query, database, params=None, compute_strategy: Literal['default'] = 'default', contact_strategy: Literal['exact', 'aprox'] = 'aprox'):
        self.query = query
        self.database = database

        self._contact_computers = {
            'default': SerialContacts
        }
        self._contact_counter = {
            'exact': ExactContacts,
            'aprox': AproxContacts
        }
        self._compute_strategy = compute_strategy
        self._contact_strategy = self._contact_counter[contact_strategy]

        self.params = params if params is not None else DEFAULT_SIM_PARAMS

    def compute(self, strategy_or_computer=None, **kwargs):
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
            contact_computer_class = self._contact_computers.get(strategy_or_computer, None)
            if contact_computer_class is None:
                strats = ', '.join(self._contact_computers.keys())
                raise ValueError(f"Unknown strategy or computer: {strategy_or_computer}. Available strategies are: {strats}.")
            contact_computer = contact_computer_class(
                self.query.universe, self.query, self.database, **kwargs
            )
        contact_computer.run(verbose=True)

        # Strategy to count and store contacts (e.g. exact, aprox, etc.)
        contact_strategy_instance = self._contact_strategy(self.query.universe, contact_computer.contact_frames, self.params.get('norm_factor'))
        contact_strategy_instance.run()

        return ComputedContacts(contact_strategy_instance, self)
