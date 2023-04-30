from typing import Callable, Literal
from prolint2.computers.contacts import ContactComputerBase, SerialContacts
from prolint2.core.typing import NestedFloatDict, NestedIterFloatDict, NestedIterIntDict

from prolint2.metrics.exact_contacts import ExactContacts
from prolint2.metrics.aprox_contacts import AproxContacts

class ContactsProvider:
    """
    Class that provides the contacts computation functionality.
    """
    def __init__(self, query, database, compute_strategy: Literal['default'] = 'default', contact_strategy: Literal['exact', 'aprox'] = 'exact'):
        self.query = query
        self.database = database

        self._contact_computers = {
            'default': SerialContacts
        }
        self._contact_counter = {
            'exact': ExactContacts,
            'aprox': AproxContacts
        }
        self._computed_contacts = None
        self._compute_strategy = compute_strategy
        self._contact_strategy = self._contact_counter[contact_strategy]

    def compute(self, strategy_or_computer=None, **kwargs):
        """
        Compute the contacts using a cythonized version of a cell-list algorithm.

        Parameters
        ----------
        strategy_or_computer : str or :class:`ContactComputerBase` ('default')
            Strategy or computer to use to compute the contacts. If a string is passed, it has to be a key in the
            **_contact_computers** dictionary. Note that only the 'default' strategy is currently available.
        kwargs : dict
            Keyword arguments to be passed to the **ContactComputerBase** class.
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
        # self._computed_contacts = contact_computer.contacts

        # Strategy to count and store contacts (e.g. exact, aprox, etc.)
        contact_strategy_instance = self._contact_strategy(self.query.universe, contact_computer.contact_frames)
        contact_strategy_instance.run()
        self._computed_contacts = contact_strategy_instance

        return self
    
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
        NestedFloatDict
            A dictionary of computed metrics for all lipids.
        """
        return self._computed_contacts.compute(metric, target_lipid_name=target_lipid_name)
    
    def apply_function(self, func: Callable, target_lipid_name=None) -> NestedFloatDict:
        """Apply a function to all lipids or a specific lipid. 

        Parameters
        ----------
        func : Callable
            The function to apply to the lipid contact durations.
        target_lipid_name : str, optional
            The name of the lipid to apply the function to. If None, the function will be applied to all lipids.

        Returns
        -------
        NestedFloatDict
            A dictionary of computed metrics for all lipids.
        """
        return self._computed_contacts.apply_function(func, target_lipid_name=target_lipid_name)
    
    @property
    def contacts(self) -> NestedIterFloatDict:
        """The computed contacts."""
        return self._computed_contacts.contacts
    
    @property
    def pooled_contacts(self) -> NestedIterFloatDict:
        """The computed contacts."""
        return self._computed_contacts.pooled_results()

    @property
    def contact_frames(self) -> NestedIterIntDict:
        """The computed contacts."""
        return self._computed_contacts.contact_frames