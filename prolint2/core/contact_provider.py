from prolint2.computers.contacts import ContactComputerBase, SerialContacts

class ContactsProvider:
    """
    Class that provides the contacts computation functionality.
    """
    def __init__(self, query, database):
        self.query = query
        self.database = database

        self._contact_computers = {
            'default': SerialContacts
        }

    def compute(self, strategy_or_computer='default', **kwargs):
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

        if isinstance(strategy_or_computer, ContactComputerBase):
            contact_computer = strategy_or_computer
        else:
            contact_computer_class = self._contact_computers[strategy_or_computer]
            contact_computer = contact_computer_class(
                self.query.universe, self.query, self.database, **kwargs
            )
        contact_computer.run(verbose=True)
        return contact_computer.contacts
