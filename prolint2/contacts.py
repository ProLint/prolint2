r"""Contacts base classes --- :mod:`prolint2.contacts`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os
import configparser
from itertools import groupby

import numpy as np
import pandas as pd

from prolint2.computers.contacts import ContactComputerBase, SerialContacts

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config.ini"))
parameters_config = config["Parameters"]


class Contacts(object):
    """Stores information to run and analyze the distance-based contacts results
    between the :class:`.prolint2.QueryProteins` and :class:`.prolint2.MembraneDatabase` groups.

    Parameters
    ----------
    query : :class:`QueryProteins`
    database : :class:`MembraneDatabase`

    Attributes
    ----------
    query : :class:`QueryProteins`
        **Query** group to use during the calculation of the contacts.
    database : :class:`MembraneDatabase`
        **Database** group to use during the calculation of the contacts.
    contacts : Array (None)
        Numpy uni-dimensional array of shape equal to the number of frames used during the calculation of the contacts.
        Each element of the array has a Scipy matrix with the pairs (i, j) defining the contacts, where *i* is the index
        of the residue in the **query** group, and *j* is the index of the residue in the **database** group. It can be populated
        using either the **compute()** or the **load()** methods.
    counts : Pandas DataFrame (None)
        Pandas DataFrame with the counted contacts. It is populated using the **count_contacts()** method.
    """

    def __init__(self, query, database):
        self.query = query
        self.database = database

        # TODO:
        # @bis: I really don't like how we have to back reference the trajectory here
        # What's the best way here? Include trajectory as an initialization argument?
        self.n_frames = query.selected.universe.trajectory.n_frames
        self.dt = self.query.selected.universe.trajectory.dt
        self.totaltime = self.query.selected.universe.trajectory.totaltime

        self._contact_computers = {
            'default': SerialContacts
            # Other contact computation strategies here
        }


    def compute(self, strategy_or_computer='default', **kwargs):
        """
        Compute the contacts using a cythonized version of a cell-list algorithm.

        Parameters
        ----------
        strategy_or_computer : str or :class:`ContactComputerBase` ('default')
            Strategy or computer to use to compute the contacts. If a string is passed, it has to be a key in the
            **_contact_computers** dictionary.
        kwargs : dict
            Keyword arguments to be passed to the **ContactComputerBase** class.
        """

        if isinstance(strategy_or_computer, ContactComputerBase):
            contact_computer = strategy_or_computer
        else:
            contact_computer_class = self._contact_computers[strategy_or_computer]
            contact_computer = contact_computer_class(
                self.query.selected.universe, self.query.selected, self.database.selected, **kwargs
            )
        contact_computer.run(verbose=True)
        return contact_computer


    # this functions allows the definition of chunks of frames with uninterrupted interactions
    # i.e. it takes a list of frames as [9, 11, 12] and it returns [1, 2]

    # def ranges(self, lst):
    #     pos = (j - i for i, j in enumerate(lst))
    #     t = 0
    #     for i, els in groupby(pos):
    #         l = len(list(els))
    #         el = lst[t]
    #         t += l
    #         yield len(range(el, el + l))

    # def computed_contacts(self):
    #     """
    #     Returns the computed contacts.
    #     """
    #     if self.contacts is None:
    #         raise ValueError("No contacts computed. Run compute() first.")
        
    #     formatted_contact_frames = {}
    #     for residue_id, lipid_dict in self.contact_frames.items():
    #         for lipid_id, frames in lipid_dict.items():
    #             contact_array = np.zeros(self.n_frames, dtype=np.int8)
    #             contact_array[frames] = 1
    #             key = (residue_id, lipid_id)
    #             formatted_contact_frames[key] = contact_array.copy()

    #     df = pd.DataFrame.from_dict(formatted_contact_frames, orient='index')
    #     df.index = pd.MultiIndex.from_tuples(df.index, names=['residueID', 'lipidID'])
    #     df = df.sort_index(
    #         level=['residueID', 'lipidID'],
    #         ascending=[True, True]
    #     )

    #     return df

    # def contacts_to_dataframe(self):
    #     """
    #     Convert the contacts dictionary to a Pandas DataFrame.

    #     Returns
    #     -------
    #     Pandas DataFrame
    #         Pandas DataFrame with all the contacts.
    #     """
    #     if not self.contacts:
    #         raise ValueError("The contacts dictionary is empty.")
    #     else:
    #         results = []
    #         keys = self.contacts.keys()
    #         for idx, protein_resi in enumerate(keys):
    #             for lip_type in self.contacts[protein_resi].keys():
    #                 for lip_res, t_frames in self.contacts[protein_resi][
    #                     lip_type
    #                 ].items():
    #                     for fr in self.contact_frames[
    #                         "{},{}".format(protein_resi, lip_res)
    #                     ]:
    #                         results.append(
    #                             (
    #                                 "Protein1",
    #                                 protein_resi,
    #                                 self.query.selected.residues[idx].resname,
    #                                 lip_type,
    #                                 lip_res,
    #                                 fr,
    #                             )
    #                         )
    #         results_df = pd.DataFrame(
    #             results,
    #             columns=[
    #                 "Protein",
    #                 "Residue ID",
    #                 "Residue Name",
    #                 "Lipid Type",
    #                 "Lipid ID",
    #                 "Frame",
    #             ],
    #         )
    #         return results_df

    # def contacts_to_metrics(self):
    #     """
    #     Convert the contacts dictionary to a Pandas DataFrame with different metrics.

    #     Returns
    #     -------
    #     Pandas DataFrame
    #         Pandas DataFrame with different metrics for the contacts.
    #     """
    #     if not self.contacts:
    #         raise ValueError("The contacts dictionary is empty.")
    #     else:
    #         metrics = []
    #         keys = self.contacts.keys()
    #         for idx, protein_resi in enumerate(keys):
    #             for lip_type in self.contacts[protein_resi].keys():
    #                 for lip_res, t_frames in self.contacts[protein_resi][
    #                     lip_type
    #                 ].items():
    #                     # getting chunks of frames with uninterrupted interactions
    #                     key = "{},{}".format(protein_resi, lip_res)
    #                     temp = list(self.ranges(self.contact_frames[key]))

    #                     # calculating metrics
    #                     metrics.append(
    #                         (
    #                             "Protein1",
    #                             protein_resi,
    #                             self.residue_names[idx],
    #                             lip_type,
    #                             lip_res,
    #                             t_frames,
    #                             t_frames / self.n_frames,
    #                             max(temp),
    #                             np.mean(temp),
    #                         )
    #                     )
    #         metrics_df = pd.DataFrame(
    #             metrics,
    #             columns=[
    #                 "Protein",
    #                 "Residue ID",
    #                 "Residue Name",
    #                 "Lipid Type",
    #                 "Lipid ID",
    #                 "Sum of all contacts",
    #                 "Occupancy",
    #                 "Longest Duration",
    #                 "Mean Duration",
    #             ],
    #         )
    #         return metrics_df

    # def export(self, filename):
    #     """
    #     Export the contacts array to a file.

    #     Parameters
    #     ----------
    #     filename : str
    #         Name of the file to export the contacts array.
    #     """
    #     print("Exporting contacts and metrics to files...")
    #     self.contacts_to_dataframe().to_csv(filename, index=False)
    #     if not isinstance(self.metrics, pd.DataFrame):            
    #         self.contacts_to_metrics().to_csv(filename.replace(".csv", "_metrics.csv"), index=False)
    #     print("Contacts successfully exported to file '{}' and metrics to '{}'!!".format(filename, filename.replace(".csv", "_metrics.csv")))

    # def filter_by_percentile(self, percentile=0.75, metric="Sum of all contacts"):
    #     """
    #     Filter the contacts by percentile.

    #     Parameters
    #     ----------
    #     percentile : float (0.75)
    #         Percentile to be used for filtering the contacts array.
    #     """
    #     if metric not in [
    #         "Sum of all contacts",
    #         "Occupancy",
    #         "Longest Duration",
    #         "Mean Duration",
    #     ]:
    #         raise ValueError("The metric is not valid.")
    #     else:
    #         return self.metrics[
    #             self.metrics[metric] > self.metrics[metric].quantile(percentile)
    #         ]
        
    # def __str__(self):
    #     if self.contacts is None:
    #         return "<prolint2.Contacts containing 0 contacts>"
    #     else:
    #         return "<prolint2.Contacts containing {} contacts>".format(len(self.contacts))

    # def __repr__(self):
    #     if self.contacts is None:
    #         return "<prolint2.Contacts containing 0 contacts>"
    #     else:
    #         return "<prolint2.Contacts containing {} contacts>".format(len(self.contacts))
