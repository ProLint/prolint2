r"""Contacts base classes --- :mod:`ufcc.contacts`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import pickle
import numpy as np
import MDAnalysis as mda
from collections import Counter
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.analysis.base import AnalysisBase
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

        # TODO:
        # @bis: I really don't like how we have to back reference the trajectory here
        # What's the best way here? Include trajectory as an initialization argument?
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
            print(
            "The parallel routine is not included in this DEMO release of the software, as our main goal is to see how the software behaves with different types and how we can improve its capabilities. Please use the serial routine for now and let us now your feedback."
            )

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


    def server_payload(self):

        # TODO:
        # protein name is hardcoded -> read protein name(s) dynamically
        # update code to handle multiple identical proteins
        # update code to handle multiple copies of different proteins
        resnames = self.query.selected.resnames
        protein_name = "QueryGroup (Protein)"
        protein = (
            protein_name  # TODO: we'll need to update this into a list and iterate over it
        )
        lipids = list(np.unique(self.database.selected.resnames))
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
        proteins = [protein_name]
        protein_counts = {protein_name: 1}

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
