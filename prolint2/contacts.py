r"""Contacts base classes --- :mod:`prolint2.contacts`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os
import configparser
from collections import Counter, defaultdict
from itertools import groupby

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase

from prolint2.metrics.metrics import create_metric
from prolint2.metrics.registries import MetricRegistry, auto_register_metrics

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "config.ini"))
parameters_config = config["Parameters"]

def process_contact_items(contact_items):
    """
    Process the contact items to get the number of contacts per residue.
    """
    processed_items = {}
    for key, value in contact_items:
        processed_items[key] = Counter(value)
    return processed_items

def transform_contacts(contacts):
    """
    Transform the contacts to a dictionary of dictionaries.

    Parameters
    ----------
    contacts : dict
        Dictionary of contacts.

    Returns
    -------
    dict: Dictionary of dictionaries.

    """
    transformed_contacts = {}
    for key, value in contacts.items():
        transformed_contacts[key] = process_contact_items(value.items())
    return transformed_contacts


class SerialContacts(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class.
    """
    # TODO:
    # @bis: The front end has the hierarch protein -> lipids -> residue
    # The data, however, are stored protein -> residue -> lipids, leading to unnecessary
    # work later on. We should modify this, so we store data in the right
    # hierarchical structure.
    def __init__(self, universe, query, database, cutoff, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)

        self.query = query
        self.database = database
        self.cutoff = cutoff

        # We need to convert to list to allow for JSON serialization
        self.q_resids = self.query.resids.tolist()
        self.db_resids = self.database.resids.tolist()
        self.db_resnames = self.database.resnames
        self.dp_resnames_unique = np.unique(self.db_resnames)

        self.contacts = self._initialize_contacts()
        self.contact_frames = defaultdict(lambda: defaultdict(list))

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _initialize_contacts(self):
        contacts = {}
        for k in self.q_resids:
            contacts[k] = {v: [] for v in self.dp_resnames_unique}
        return contacts

    def _get_residue_lipid_info(self, pair):
        residue_id = self.q_resids[pair[0]]
        lipid_id = self.db_resids[pair[1]]
        lipid_name = self.db_resnames[pair[1]]
        return residue_id, lipid_id, lipid_name

    def _compute_pairs(self):
        """
        Compute the pairs of residues and lipids that are within the cutoff distance.
        """
        gridsearch = FastNS(
            self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True
        )
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()
        return pairs

    def _single_frame(self):
        """
        Compute the contacts for a single frame.
        """
        pairs = self._compute_pairs()

        existing_pairs = defaultdict(set)
        for pair in pairs:
            residue_lipid = self._get_residue_lipid_info(pair)
            if residue_lipid not in existing_pairs:

                existing_pairs[residue_lipid].add(residue_lipid)

                residue_id, lipid_id, lipid_name = residue_lipid
                self.contacts[residue_id][lipid_name].append(lipid_id)
                self.contact_frames[residue_id][lipid_id].append(self._frame_index)

    def _conclude(self):
        self.contacts = transform_contacts(self.contacts)

class SerialDistances(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class.
    """
    # TODO:
    # @bis: The front end has the hierarch protein -> lipids -> residue
    # The data, however, are stored protein -> residue -> lipids, leading to unnecessary
    # work later on. We should modify this, so we store data in the right
    # hierarchical structure.
    def __init__(
        self, universe, query, database, lipid_id, residue_id, frame_filter, **kwargs
    ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.frame_filter = frame_filter
        frame_range = np.arange(len(self.frame_filter))
        self.frame_mapping = {k: v for k, v in zip(self.frame_filter, frame_range)}

        self.lipid_atomgroup = self.database.select_atoms(f"resid {lipid_id}")
        self.resid_atomgroup = self.query.select_atoms(f"resid {residue_id}")
        self.lipid_atomnames = self.lipid_atomgroup.names.tolist()
        self.resid_atomnames = self.resid_atomgroup.names.tolist()
        self.result_array = None
        self.distance_array = None

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

    def _prepare(self):
        self.result_array = np.zeros(
            (
                len(self.frame_filter),
                self.lipid_atomgroup.n_atoms,
                self.resid_atomgroup.n_atoms,
            )
        )

    def _single_frame(self):
        if self._frame_index in self.frame_filter:
            r = distances.distance_array(
                self.lipid_atomgroup.positions,
                self.resid_atomgroup.positions,
                box=self.database.universe.dimensions,
            )
            # print ('frame iterator: ', self._frame_index)
            self.result_array[self.frame_mapping[self._frame_index]] = r

    def _conclude(self):
        self.distance_array = np.mean(self.result_array, axis=0)
        del self.result_array


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
        self.residue_names = self.query.selected.residues.resnames
        self.residue_ids = self.query.selected.residues.resids
        self.cutoff = None
        self.contacts = None
        self.contact_frames = None
        self.metrics = None

        self.registry = MetricRegistry()
        self._register_metrics()

        # TODO:
        # @bis: I really don't like how we have to back reference the trajectory here
        # What's the best way here? Include trajectory as an initialization argument?
        self.n_frames = query.selected.universe.trajectory.n_frames
        self.dt = self.query.selected.universe.trajectory.dt
        self.totaltime = self.query.selected.universe.trajectory.totaltime

    def _register_metrics(self):
        auto_register_metrics(self.registry, 'prolint2.metrics.metrics')

    def compute(self, cutoff=int(parameters_config["cutoff"]), get_metrics=False):
        """
        Compute the cutoff distance-based contacts using a cythonized version of a cell-list algorithm.

        Parameters
        ----------
        cutoff : int (7)
            Value in Angstrom to be used as cutoff for the calculation of the contacts.
        """
        self.cutoff = cutoff
        assert isinstance(
            self.query.selected,
            (mda.core.groups.AtomGroup),
        ), "the query has to be an AtomGroup"
        assert isinstance(
            self.database.selected,
            (mda.core.groups.AtomGroup),
        ), "the database has to be an AtomGroup"
        temp_instance = SerialContacts(
            self.query.selected.universe,
            self.query.selected,
            self.database.selected,
            cutoff,
        )
        temp_instance.run(verbose=True)

        self.contacts = temp_instance.contacts
        self.contact_frames = temp_instance.contact_frames
        if get_metrics:
            self.metrics = self.contacts_to_metrics()

    # this functions allows the definition of chunks of frames with uninterrupted interactions
    # i.e. it takes a list of frames as [9, 11, 12] and it returns [1, 2]
    def ranges(self, lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield len(range(el, el + l))

    def computed_contacts(self):
        """
        Returns the computed contacts.
        """
        if self.contacts is None:
            raise ValueError("No contacts computed. Run compute() first.")
        
        formatted_contact_frames = {}
        for residue_id, lipid_dict in self.contact_frames.items():
            for lipid_id, frames in lipid_dict.items():
                contact_array = np.zeros(self.n_frames, dtype=np.int8)
                contact_array[frames] = 1
                key = (residue_id, lipid_id)
                formatted_contact_frames[key] = contact_array.copy()

        df = pd.DataFrame.from_dict(formatted_contact_frames, orient='index')
        df.index = pd.MultiIndex.from_tuples(df.index, names=['residueID', 'lipidID'])
        df = df.sort_index(
            level=['residueID', 'lipidID'],
            ascending=[True, True]
        )

        return df

    def contacts_to_dataframe(self):
        """
        Convert the contacts dictionary to a Pandas DataFrame.

        Returns
        -------
        Pandas DataFrame
            Pandas DataFrame with all the contacts.
        """
        if not self.contacts:
            raise ValueError("The contacts dictionary is empty.")
        else:
            results = []
            keys = self.contacts.keys()
            for idx, protein_resi in enumerate(keys):
                for lip_type in self.contacts[protein_resi].keys():
                    for lip_res, t_frames in self.contacts[protein_resi][
                        lip_type
                    ].items():
                        for fr in self.contact_frames[
                            "{},{}".format(protein_resi, lip_res)
                        ]:
                            results.append(
                                (
                                    "Protein1",
                                    protein_resi,
                                    self.query.selected.residues[idx].resname,
                                    lip_type,
                                    lip_res,
                                    fr,
                                )
                            )
            results_df = pd.DataFrame(
                results,
                columns=[
                    "Protein",
                    "Residue ID",
                    "Residue Name",
                    "Lipid Type",
                    "Lipid ID",
                    "Frame",
                ],
            )
            return results_df

    def contacts_to_metrics(self):
        """
        Convert the contacts dictionary to a Pandas DataFrame with different metrics.

        Returns
        -------
        Pandas DataFrame
            Pandas DataFrame with different metrics for the contacts.
        """
        if not self.contacts:
            raise ValueError("The contacts dictionary is empty.")
        else:
            metrics = []
            keys = self.contacts.keys()
            for idx, protein_resi in enumerate(keys):
                for lip_type in self.contacts[protein_resi].keys():
                    for lip_res, t_frames in self.contacts[protein_resi][
                        lip_type
                    ].items():
                        # getting chunks of frames with uninterrupted interactions
                        key = "{},{}".format(protein_resi, lip_res)
                        temp = list(self.ranges(self.contact_frames[key]))

                        # calculating metrics
                        metrics.append(
                            (
                                "Protein1",
                                protein_resi,
                                self.residue_names[idx],
                                lip_type,
                                lip_res,
                                t_frames,
                                t_frames / self.n_frames,
                                max(temp),
                                np.mean(temp),
                            )
                        )
            metrics_df = pd.DataFrame(
                metrics,
                columns=[
                    "Protein",
                    "Residue ID",
                    "Residue Name",
                    "Lipid Type",
                    "Lipid ID",
                    "Sum of all contacts",
                    "Occupancy",
                    "Longest Duration",
                    "Mean Duration",
                ],
            )
            return metrics_df

    def export(self, filename):
        """
        Export the contacts array to a file.

        Parameters
        ----------
        filename : str
            Name of the file to export the contacts array.
        """
        print("Exporting contacts and metrics to files...")
        self.contacts_to_dataframe().to_csv(filename, index=False)
        if not isinstance(self.metrics, pd.DataFrame):            
            self.contacts_to_metrics().to_csv(filename.replace(".csv", "_metrics.csv"), index=False)
        print("Contacts successfully exported to file '{}' and metrics to '{}'!!".format(filename, filename.replace(".csv", "_metrics.csv")))

    def filter_by_percentile(self, percentile=0.75, metric="Sum of all contacts"):
        """
        Filter the contacts by percentile.

        Parameters
        ----------
        percentile : float (0.75)
            Percentile to be used for filtering the contacts array.
        """
        if metric not in [
            "Sum of all contacts",
            "Occupancy",
            "Longest Duration",
            "Mean Duration",
        ]:
            raise ValueError("The metric is not valid.")
        else:
            return self.metrics[
                self.metrics[metric] > self.metrics[metric].quantile(percentile)
            ]
        
    def server_payload(self, metric="max", custom_user_function=None):

        # TODO:
        # protein name is hardcoded -> read protein name(s) dynamically
        # update code to handle multiple identical proteins
        # update code to handle multiple copies of different proteins
        protein_name = "Protein" # TODO: we'll need to update this into a list and iterate over it
        proteins = [protein_name]
        protein_counts = {protein_name: 1}

        # residue_contacts = {}
        metric_instance = create_metric(
            self.contacts, 
            metrics=[metric], 
            metric_registry=self.registry, 
            output_format='dashboard', 
            lipid_type=self.database.lipid_types().tolist()[0], 
            residue_names=self.residue_names, 
            residue_ids=self.residue_ids
        )
        residue_contacts = metric_instance.compute(dt=self.dt, totaltime=self.totaltime)

        lipid_counts = self.database.lipid_count()
        total_lipid_sum = sum(lipid_counts.values())
        sub_data = []
        for lipid, count in lipid_counts.items():
            sub_data.append({"category": lipid, "value": "{:.2f}".format(count / total_lipid_sum)})

        pie_data = []
        for protein in proteins:
            value = protein_counts[protein] / sum(protein_counts.values())

            protein_pdata = {
                "category": protein_name,
                "value": "{:.2f}".format(value),
                "subData": sub_data,
            }
            pie_data.append(protein_pdata)

        payload = {
            "data": {protein_name: residue_contacts},
            "proteins": [protein_name],
            "lipids": self.database.lipid_types().tolist(),
            "pie_data": pie_data,  # TODO: include protein info
        }

        return payload

    def __str__(self):
        if self.contacts is None:
            return "<prolint2.Contacts containing 0 contacts>"
        else:
            return "<prolint2.Contacts containing {} contacts>".format(len(self.contacts))

    def __repr__(self):
        if self.contacts is None:
            return "<prolint2.Contacts containing 0 contacts>"
        else:
            return "<prolint2.Contacts containing {} contacts>".format(len(self.contacts))
