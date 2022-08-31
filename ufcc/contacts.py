r"""Contacts base classes --- :mod:`ufcc.contacts`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import pickle
import numpy as np
import pandas as pd
import scipy.stats
import scipy.sparse
from tqdm import tqdm
import MDAnalysis as mda
from collections import Counter
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.analysis.base import AnalysisBase
from .parallel import ParallelAnalysisBase

from MDAnalysis.analysis import distances


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
        self.q_resids = self.query.resindices.tolist()
        self.db_resids = self.database.resindices.tolist()
        self.db_resnames = self.database.resnames
        self.dp_resnames_unique = np.unique(self.db_resnames)

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _prepare(self):
        print("PREPARING contact_frames")
        self.contacts = {
            k: {v: [] for v in self.dp_resnames_unique} for k in self.q_resids
        }
        self.contacts_sum = {
            k: {v: 0 for v in self.dp_resnames_unique} for k in self.q_resids
        }
        self.contact_frames = {}

    def _single_frame(self):
        gridsearch = FastNS(
            self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True
        )
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        existing_pairs = {}
        for p in pairs:
            residue_id = self.q_resids[p[0]]
            lipid_id = self.db_resids[p[1]]
            string = f"{residue_id},{lipid_id}"

            # if self._frame_index == 0 and residue_id < 200:
            #     print (p[0], p[1], residue_id, lipid_id)
            # NOTE:
            # We want to keep track of frames the cutoff is satisfied
            # and also the pairs that satisfied the cutoff -> this can be used to avoid
            # the distance array analysis necessary later.
            # frame_pairs = (self._frame_index, p)
            # if string in self.contact_frames:
            #     self.contact_frames[string].append(frame_pairs)
            # else:
            #     self.contact_frames[string] = [frame_pairs]

            if f"{residue_id}{lipid_id}" in existing_pairs:
                continue
            existing_pairs[f"{residue_id}{lipid_id}"] = True

            # TODO:
            # @bis: we may be able to get further performance improvements by
            # using the Counter object with its update methods.

            # TODO:
            # these IDs are not guaranteed to be unique:
            # For systems containing multiple proteins
            # For very large systems with duplicate lipid residue IDs (e.g. two instances of 1234CHOL)
            lipid_name = self.db_resnames[p[1]]
            self.contacts_sum[residue_id][lipid_name] += 1
            self.contacts[residue_id][lipid_name].append(lipid_id)

            # NOTE:
            # We want to keep track of frames the cutoff is satisfied
            # the self.contact_frames dict gets very large and may not be feasible for large systems.
            # In general, it's not a method that's going to scale well. Given the backend we have, it
            # makes more sense to store results in a temporary SQL database. Retrieval will be superfast,
            # and we can do much more that way.
            if string in self.contact_frames:
                self.contact_frames[string].append(self._frame_index)
            else:
                self.contact_frames[string] = [self._frame_index]

    def _conclude(self):
        # self.contacts_sum = dict(map(lambda x: (x[0], Counter(x[1])), self.contacts_sum.items()))
        self.contacts = dict(
            map(
                lambda x: (
                    x[0],
                    dict(map(lambda y: (y[0], Counter(y[1])), x[1].items())),
                ),
                self.contacts.items(),
            )
        )


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
        self.lipid_atomgroup = self.database.select_atoms(f"resid {lipid_id+1}")
        self.resid_atomgroup = self.query.select_atoms(f"resid {residue_id+1}")
        self.lipid_atomnames = self.lipid_atomgroup.names.tolist()
        self.resid_atomnames = self.resid_atomgroup.names.tolist()

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

    def _prepare(self):
        self.result_array = np.zeros(
            (self.n_frames, self.lipid_atomgroup.n_atoms, self.resid_atomgroup.n_atoms)
        )

    def _single_frame(self):
        if self._frame_index in self.frame_filter:
            r = distances.distance_array(
                self.lipid_atomgroup.positions,
                self.resid_atomgroup.positions,
                box=self.database.universe.dimensions,
            )
            self.result_array[self._frame_index] = r

    def _conclude(self):
        self.distance_array = np.mean(self.result_array, axis=0)
        del self.result_array


class ParallelContacts(ParallelAnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *parallel* approach.
    """

    def __init__(self, universe, query, database, cutoff, **kwargs):

        super().__init__(universe.universe.trajectory, (query, database))
        self.query = query
        self.query_n_residues = self.query.n_residues
        self.database = database
        self.database_n_residues = self.database.n_residues
        self.cutoff = cutoff

        # to allow for non-sequential resindices
        self._sorted_protein_resindices = (
            scipy.stats.rankdata(self.query.resindices, method="dense") - 1
        )
        self._sorted_membrane_resindices = (
            scipy.stats.rankdata(self.database.resindices, method="dense") - 1
        )

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _prepare(self):
        # Initialize empty np.array with results
        self.contacts = None

    def _single_frame(self, ts, atomgroups):
        # Get the results and populate the results dictionary
        gridsearch = FastNS(
            self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True
        )
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        query_residx, database_residx = np.unique(
            np.array(
                [
                    [
                        self._sorted_protein_resindices[pair[0]],
                        self._sorted_membrane_resindices[pair[1]],
                    ]
                    for pair in pairs
                ]
            ),
            axis=0,
        ).T

        # store neighbours for this frame
        data = np.ones_like(query_residx)
        return (
            ts.frame,
            scipy.sparse.csr_matrix(
                (data, (query_residx, database_residx)),
                dtype=np.int8,
                shape=(self.query_n_residues, self.database_n_residues),
            ),
        )

    def _conclude(self):
        self.contacts = np.array(
            [l[1] for l in sorted(np.vstack(self._results), key=lambda tup: tup[0])]
        )


class Runner(object):
    """
    Class to configure the runner for the calculations of distances. As the *parallel* routine uses
    the parallel computing library **Dask**, that can be setted up to run on remotes machines. The aim of
    this runner class is to define the variables needed for the Dask scheduler, but so far this a very
    simple class that has the attributes below to run on a single local machine:

    Attributes
    ----------
    backend : str (*serial*)
        Backend to run the contacts calculation (can be either *serial* or *parallel*).
    n_jobs : int (-1)
        Number of cores to use with the *parallel* backend. By default **ufcc** will use all the cores.
    """

    def __init__(self):
        self.backend = "serial"
        self.n_jobs = -1
        # TODO
        # add funcionalities to run analysis on HPC machines


class Contacts(object):
    """Stores information to run and analyze the distance-based contacts results
    between the :class:`.ufcc.QueryProteins` and :class:`.ufcc.MembraneDatabase` groups.

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
    runner : :class:`Runner`
        Runner object to define the backend (*serial* or *parallel*) and n_jobs to use during the calculation of the contacts.
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
        self.runner = Runner()
        self.cutoff = None
        self.contacts = None
        self.counts = None
        self.contact_metrics = None

        # TODO:
        # @bis: I really don't like how we have to back reference the trajectory here
        # What's the best way here? Include trajectory as an initialization argument?
        self.database_unique = np.unique(self.database.selected.resnames)
        self.dt = self.query.selected.universe.trajectory.dt
        self.totaltime = self.query.selected.universe.trajectory.totaltime

    def compute(self, cutoff=7):
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
        # TODO:
        # @bis: store backend list in the project config file.
        if self.runner.backend == None or self.runner.backend not in [
            "serial",
            "parallel",
        ]:
            raise ValueError(
                "You have to select a proper backend before running the contacts routine. \n Valid options: 'serial', 'parallel'"
            )
        if self.runner.backend == "serial":
            temp_instance = SerialContacts(
                self.query.selected.universe,
                self.query.selected,
                self.database.selected,
                cutoff,
            )
            temp_instance.run(verbose=True)
        elif self.runner.backend == "parallel":
            temp_instance = ParallelContacts(
                self.query.selected.universe,
                self.query.selected,
                self.database.selected,
                cutoff,
            )
            temp_instance.run(n_jobs=self.runner.n_jobs)

        self.contacts = temp_instance.contacts
        self.contacts_sum = temp_instance.contacts_sum
        self.contact_frames = temp_instance.contact_frames

    def save(self, path="contacts.pkl"):
        """
        Store the contacts information in a pickle file for later usage.

        Parameters
        ----------
        path : file
            Path to file to save the contacts information. ('contacts.pkl')
        """
        with open(path, "wb") as f:
            pickle.dump(self.contacts, f)

    def load(self, path="contacts.pkl"):
        """
        Load the contacts information from a pickle file.

        Parameters
        ----------
        path : file
            Path to file to load the contacts information from. ('contacts.pkl')
        """
        with open(path, "rb") as f:
            self.contacts = pickle.load(f)

    def count_contacts(self):
        """
        Count the number of each contact type at each frame.
        """

        if self.contacts is None:
            raise ValueError(
                "contacts attribute is None: use .compute() before calling .count_neighbours()"
            )

        # Use lipid resnames to distinguish lipids
        count_by = np.full(
            (
                self.database.selected.n_residues,
                self.database.selected.universe.trajectory.n_frames,
            ),
            fill_value=self.database.selected.residues.resnames[:, np.newaxis],
        )
        count_by_labels = {
            label: index
            for index, label in enumerate(
                np.unique(self.database.selected.residues.resnames)
            )
        }

        # create output array
        all_counts = np.full(
            (
                self.query.selected.n_residues,
                self.query.selected.universe.trajectory.n_frames,
                len(count_by_labels),
            ),
            fill_value=0,
            dtype=np.uint8,  # count can't be negative, and no lipid will have more than 255 neighbours
        )

        # For counts we need to know which column of the output array to add counts to for each lipid type
        type_index = {value: index for index, value in enumerate(count_by_labels)}

        # Get counts at each frame
        for frame_index, contacts in tqdm(
            enumerate(self.contacts),
            total=self.query.selected.universe.trajectory.n_frames,
        ):
            ref, neigh = contacts.nonzero()
            unique, counts = np.unique(
                [ref, [type_index[t] for t in count_by[neigh, frame_index]]],
                axis=1,
                return_counts=True,
            )

            r, t = unique  # reference index (r) and type index (t)
            all_counts[r, frame_index, t] = counts

        labels = np.full(
            (
                self.query.selected.n_residues,
                self.query.selected.universe.trajectory.n_frames,
            ),
            fill_value=self.query.selected.residues.macros[:, np.newaxis],
        )
        labels = labels.reshape(
            self.query.selected.n_residues
            * self.query.selected.universe.trajectory.n_frames
        )

        # Assemble data for the DataFrame
        residue_labels = np.full(
            (
                self.query.selected.n_residues,
                self.query.selected.universe.trajectory.n_frames,
            ),
            fill_value=self.query.selected.residues.resnames[:, np.newaxis],
        )
        residue_labels = residue_labels.reshape(
            self.query.selected.n_residues
            * self.query.selected.universe.trajectory.n_frames
        )
        # labels = np.array([list(count_by_labels)[type_index[frame_index]] for lipid in count_by for frame_index in lipid])

        resindices = np.full(
            (
                self.query.selected.n_residues,
                self.query.selected.universe.trajectory.n_frames,
            ),
            fill_value=self.query.selected.residues.resindices[:, np.newaxis],
        )
        resindices = resindices.reshape(
            self.query.selected.n_residues
            * self.query.selected.universe.trajectory.n_frames
        )

        frames = np.full(
            (
                self.query.selected.n_residues,
                self.query.selected.universe.trajectory.n_frames,
            ),
            fill_value=range(self.query.selected.universe.trajectory.n_frames),
        )
        frames = frames.reshape(
            self.query.selected.n_residues
            * self.query.selected.universe.trajectory.n_frames
        )

        all_counts = all_counts.reshape(
            self.query.selected.n_residues
            * self.query.selected.universe.trajectory.n_frames,
            len(count_by_labels),
        )
        total_counts = np.sum(all_counts, axis=1)

        # Create the dataframe
        # counts = pd.DataFrame(data=residue_labels, columns=["Residue"])
        counts = pd.DataFrame(data=labels, columns=["Protein"])

        counts["Residue"] = residue_labels
        counts["ResID"] = resindices
        counts["FrameID"] = frames

        for count_by_label in count_by_labels:
            counts[f"# {count_by_label}"] = all_counts.T[type_index[count_by_label]]

        counts["Total"] = total_counts

        # make every column except the label take on integer values
        for column in counts.columns[2:]:
            counts[column] = pd.to_numeric(counts[column])

        self.counts = counts

    def get_metrics(self, save_file="", server=False):
        if self.contacts is None or self.counts is None:
            raise ValueError(
                "contacts or counts attributes are None: use .compute() and .count_neighbours() methods  before calling get_metrics"
            )

        n_frames = self.database.selected.universe.trajectory.n_frames
        unique_lip_resnames = np.unique(self.database.selected.residues.resnames)
        protein = list(self.counts[self.counts["FrameID"] == 0]["Protein"]) * len(
            unique_lip_resnames
        )
        res_ids = list(self.counts[self.counts["FrameID"] == 0]["ResID"]) * len(
            unique_lip_resnames
        )
        res_names = list(self.counts[self.counts["FrameID"] == 0]["Residue"]) * len(
            unique_lip_resnames
        )
        radius = self.cutoff
        lipids = []
        occupancy = []
        lipid_number = []
        sum_of_all_contacts = []
        for lip in unique_lip_resnames:
            for res in list(self.counts[self.counts["FrameID"] == 0]["ResID"]):
                lipids.append(lip)
            pivot_t = self.counts.pivot_table(
                index="ResID", columns="FrameID", values=f"# {lip}"
            )
            occupancy += list(pivot_t.astype(bool).sum(axis=1) * 100 / n_frames)
            lipid_number += list(
                pivot_t.sum(axis=1)
                / np.count_nonzero(self.database.selected.residues.resnames == lip)
            )
            # sum_of_all_contacts += list(pivot_t.sum(axis=1) / n_frames)

        contact_metrics = pd.DataFrame({"Protein": protein})
        contact_metrics["ResID"] = res_ids
        contact_metrics["ResName"] = res_names
        contact_metrics["Lipids"] = lipids
        contact_metrics["Radius"] = radius
        contact_metrics["Occupancy"] = occupancy
        contact_metrics["Lipid_Number"] = lipid_number
        # contact_metrics['Sum_of_all_Contacts'] = sum_of_all_contacts

        if save_file != "":
            contact_metrics.to_csv(save_file, index=False)

        self.contact_metrics = contact_metrics

        if server:
            return contact_metrics

    def server_payload(self):

        # TODO:
        # protein name is hardcoded -> read protein name(s) dynamically
        # update code to handle multiple identical proteins
        # update code to handle multiple copies of different proteins
        resnames = self.query.selected.resnames
        protein = (
            "GIRK"  # TODO: we'll need to update this into a list and iterate over it
        )
        lipids = list(self.database_unique)
        sub_data = {
            k: {"category": k, "value": 0} for k in lipids
        }  # TODO: we need to generate sub_data for each protein.
        js = {protein: {k: [] for k in lipids}}
        for residue, contact_counter in self.contacts_sum.items():
            for lipid, contact_sum in contact_counter.items():
                sub_data[lipid]["value"] += contact_sum
                metric = (
                    contact_sum * self.dt
                ) / self.totaltime  # TODO: do we have to substract 1 frame here?
                js[protein][lipid].append(
                    {
                        "residue": f"{resnames[residue]} {residue+1}",
                        "value": float("{:.2f}".format(metric)),
                    }
                )

        sub_data = list(sub_data.values())
        norm_with = sum([x["value"] for x in sub_data])
        sub_data = [
            {
                "category": d["category"],
                "value": "{:.2f}".format(d["value"] / norm_with),
            }
            for d in sub_data
        ]

        # return js, {protein: sub_data}

        # TODO:
        # Hardcoded
        proteins = ["GIRK"]
        protein_counts = {"GIRK": 1}

        pie_data = []
        for protein in proteins:
            value = protein_counts[protein] / sum(protein_counts.values())

            protein_pdata = {
                "category": protein,
                "value": "{:.2f}".format(value),
                "subData": sub_data,
            }
            pie_data.append(protein_pdata)

        # ganttApp toy data
        gantt_data = [
            {
                "category": "Lipid 1",
                "startFrame": 0,
                "endFrame": 10,
            },
            {
                "category": "Lipid 1",
                "startFrame": 45,
                "endFrame": 75,
            },
            {
                "category": "Lipid 1",
                "startFrame": 90,
                "endFrame": 100,
            },
            {
                "category": "Lipid 2",
                "startFrame": 10,
                "endFrame": 35,
            },
            {
                "category": "Lipid 2",
                "startFrame": 30,
                "endFrame": 60,
            },
        ]
        top_10_lipids = ["Lipid 1", "Lipid 2"]

        # payload should include the entire data. The backend can process it then based on client requests
        payload = {
            "data": js,
            "proteins": [protein],
            "lipids": lipids,
            "pie_data": pie_data,  # TODO: include protein info
            "gantt_data": gantt_data,
            "top_10_lipids": top_10_lipids,
        }

        return payload

    def __str__(self):
        # TODO:
        # @bis: This is not working with the 'prolint_serial' backend, bc
        # it is not an isntance of ndarray.
        if not isinstance(self.contacts, np.ndarray):
            # TODO: should this say that contacts needs to be initialized?
            # one edge case here would be cases when there really are no contacts
            return "<ufcc.Contacts containing 0 contacts>"
        else:
            return "<ufcc.Contacts containing {} contacts>".format(len(self.contacts))

    def __repr__(self):
        # TODO:
        # @bis: same as above.
        if not isinstance(self.contacts, np.ndarray):
            return "<ufcc.Contacts containing 0 contacts>"
        else:
            return "<ufcc.Contacts containing {} contacts>".format(len(self.contacts))
