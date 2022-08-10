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
from collections import defaultdict, Counter
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.analysis.base import AnalysisBase
from .parallel import ParallelAnalysisBase
from. w2plp import LPContacts


class SerialContacts(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It heritages from the MDAnalysis AnalysisBase class.
    """

    def __init__(self, universe, query, database, cutoff, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.cutoff = cutoff

        # to allow for non-sequential resindices
        self._sorted_protein_resindices = scipy.stats.rankdata(self.query.resindices, method="dense") - 1
        self._sorted_membrane_resindices = scipy.stats.rankdata(self.database.resindices, method="dense") - 1

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _prepare(self):
        # Initialize empty np.array with results
        self.contacts = np.zeros(self.n_frames, dtype=object)

    def _single_frame(self):
        # Get the results and populate the results dictionary
        gridsearch = FastNS(self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True)
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        query_residx, database_residx = np.unique(np.array(
            [[self._sorted_protein_resindices[pair[0]], self._sorted_membrane_resindices[pair[1]]] for pair in pairs]),
                                                  axis=0).T

        # store neighbours for this frame
        data = np.ones_like(query_residx)
        self.contacts[self._frame_index] = scipy.sparse.csr_matrix(
            (data, (query_residx, database_residx)),
            dtype=np.int8,
            shape=(self.query.n_residues, self.database.n_residues))

    # def _conclude(self):
    #     # OPTIONAL
    #     # Called once iteration on the trajectory is finished.
    #     # Apply normalisation and averaging to results here.
    #     self.result = np.asarray(self.result) / np.sum(self.result)

class ProLintSerialContacts(AnalysisBase):
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
        self.q_resids = self.query.resindices
        self.db_resids = self.database.resindices
        self.db_resnames = self.database.resnames
        self.dp_resnames_unique = np.unique(self.db_resnames)

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _prepare(self):
        self.contacts = {k: {v:[] for v in self.dp_resnames_unique} for k in self.q_resids}
        self.contacts_sum = {k: {v: 0 for v in self.dp_resnames_unique} for k in self.q_resids}

    def _single_frame(self):
        gridsearch = FastNS(self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True)
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        existing_pairs = {}
        for p in pairs:
            residue_id = self.q_resids[p[0]]
            lipid_id = self.db_resids[p[1]]

            if f'{residue_id}{lipid_id}' in existing_pairs: continue
            existing_pairs[f'{residue_id}{lipid_id}'] = True

            # TODO:
            # @bis: we may be able to get further performance improvements by
            # using the Counter object with its update methods.
            lipid_name = self.db_resnames[p[1]]
            self.contacts_sum[residue_id][lipid_name] += 1
            self.contacts[residue_id][lipid_name].append(lipid_id)

    def _conclude(self):
        # self.contacts_sum = dict(map(lambda x: (x[0], Counter(x[1])), self.contacts_sum.items()))
        self.contacts = dict(map(
            lambda x: (x[0], dict(map(
                lambda y: (y[0], Counter(y[1])), x[1].items())
                )), self.contacts.items()))


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
        self._sorted_protein_resindices = scipy.stats.rankdata(self.query.resindices, method="dense") - 1
        self._sorted_membrane_resindices = scipy.stats.rankdata(self.database.resindices, method="dense") - 1

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
        gridsearch = FastNS(self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True)
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        query_residx, database_residx = np.unique(np.array(
            [[self._sorted_protein_resindices[pair[0]], self._sorted_membrane_resindices[pair[1]]] for pair in pairs]),
                                                  axis=0).T

        # store neighbours for this frame
        data = np.ones_like(query_residx)
        return (ts.frame,
                scipy.sparse.csr_matrix((data, (query_residx, database_residx)),
                                        dtype=np.int8,
                                        shape=(self.query_n_residues, self.database_n_residues)))

    def _conclude(self):
        self.contacts = np.array([l[1] for l in sorted(np.vstack(self._results), key=lambda tup: tup[0])])


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
        self.backend = 'serial'
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
        if self.runner.backend == None or self.runner.backend not in ['serial', 'parallel', 'prolint_serial']:
            raise ValueError(
                "You have to select a proper backend before running the contacts routine. \n Valid options: 'serial', 'parallel'"
            )
        if self.runner.backend == 'serial':
            temp_instance = SerialContacts(self.query.selected.universe, self.query.selected, self.database.selected,
                                           cutoff)
            temp_instance.run(verbose=True)
        elif self.runner.backend == 'prolint_serial':
            temp_instance = ProLintSerialContacts(self.query.selected.universe, self.query.selected, self.database.selected, cutoff)
            temp_instance.run(verbose=True)
        elif self.runner.backend == 'parallel':
            temp_instance = ParallelContacts(self.query.selected.universe, self.query.selected, self.database.selected,
                                             cutoff)
            temp_instance.run(n_jobs=self.runner.n_jobs)

        self.contacts = temp_instance.contacts
        self.contacts_sum = temp_instance.contacts_sum

    def save(self, path='contacts.pkl'):
        """
        Store the contacts information in a pickle file for later usage.

        Parameters
        ----------
        path : file
            Path to file to save the contacts information. ('contacts.pkl')
        """
        with open(path, 'wb') as f:
            pickle.dump(self.contacts, f)

    def load(self, path='contacts.pkl'):
        """
        Load the contacts information from a pickle file.

        Parameters
        ----------
        path : file
            Path to file to load the contacts information from. ('contacts.pkl')
        """
        with open(path, 'rb') as f:
            self.contacts = pickle.load(f)

    def count_contacts(self):
        """
        Count the number of each contact type at each frame.
        """

        if self.contacts is None:
            raise ValueError("contacts attribute is None: use .compute() before calling .count_neighbours()")

        # Use lipid resnames to distinguish lipids
        count_by = np.full(
            (self.database.selected.n_residues, self.database.selected.universe.trajectory.n_frames),
            fill_value=self.database.selected.residues.resnames[:, np.newaxis],
        )
        count_by_labels = {
            label: index
            for index, label in enumerate(np.unique(self.database.selected.residues.resnames))
        }

        # create output array
        all_counts = np.full(
            (self.query.selected.n_residues, self.query.selected.universe.trajectory.n_frames, len(count_by_labels)),
            fill_value=0,
            dtype=np.uint8  # count can't be negative, and no lipid will have more than 255 neighbours
        )

        # For counts we need to know which column of the output array to add counts to for each lipid type
        type_index = {value: index for index, value in enumerate(count_by_labels)}

        # Get counts at each frame
        for frame_index, contacts in tqdm(enumerate(self.contacts),
                                          total=self.query.selected.universe.trajectory.n_frames):
            ref, neigh = contacts.nonzero()
            unique, counts = np.unique([ref, [type_index[t] for t in count_by[neigh, frame_index]]],
                                       axis=1,
                                       return_counts=True)

            r, t = unique  # reference index (r) and type index (t)
            all_counts[r, frame_index, t] = counts

        labels = np.full((self.query.selected.n_residues, self.query.selected.universe.trajectory.n_frames),
                         fill_value=self.query.selected.residues.macros[:, np.newaxis])
        labels = labels.reshape(self.query.selected.n_residues * self.query.selected.universe.trajectory.n_frames)

        # Assemble data for the DataFrame
        residue_labels = np.full((self.query.selected.n_residues, self.query.selected.universe.trajectory.n_frames),
                                 fill_value=self.query.selected.residues.resnames[:, np.newaxis])
        residue_labels = residue_labels.reshape(self.query.selected.n_residues *
                                                self.query.selected.universe.trajectory.n_frames)
        # labels = np.array([list(count_by_labels)[type_index[frame_index]] for lipid in count_by for frame_index in lipid])

        resindices = np.full((self.query.selected.n_residues, self.query.selected.universe.trajectory.n_frames),
                             fill_value=self.query.selected.residues.resindices[:, np.newaxis])
        resindices = resindices.reshape(self.query.selected.n_residues *
                                        self.query.selected.universe.trajectory.n_frames)

        frames = np.full((self.query.selected.n_residues, self.query.selected.universe.trajectory.n_frames),
                         fill_value=range(self.query.selected.universe.trajectory.n_frames))
        frames = frames.reshape(self.query.selected.n_residues * self.query.selected.universe.trajectory.n_frames)

        all_counts = all_counts.reshape(
            self.query.selected.n_residues * self.query.selected.universe.trajectory.n_frames, len(count_by_labels))
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

    def get_metrics(self, save_file='', server=False):
        if self.contacts is None or self.counts is None:
            raise ValueError("contacts or counts attributes are None: use .compute() and .count_neighbours() methods  before calling get_metrics")

        n_frames = self.database.selected.universe.trajectory.n_frames
        unique_lip_resnames = np.unique(self.database.selected.residues.resnames)
        protein = list(self.counts[self.counts['FrameID'] == 0]['Protein']) * len(unique_lip_resnames)
        res_ids = list(self.counts[self.counts['FrameID'] == 0]['ResID']) * len(unique_lip_resnames)
        res_names = list(self.counts[self.counts['FrameID'] == 0]['Residue']) * len(unique_lip_resnames)
        radius = self.cutoff
        lipids = []
        occupancy = []
        lipid_number = []
        sum_of_all_contacts = []
        for lip in unique_lip_resnames:
            for res in list(self.counts[self.counts['FrameID'] == 0]['ResID']):
                lipids.append(lip)
            pivot_t = self.counts.pivot_table(index='ResID', columns='FrameID', values=f"# {lip}")
            occupancy += list(pivot_t.astype(bool).sum(axis=1) * 100 / n_frames)
            lipid_number += list(pivot_t.sum(axis=1) / np.count_nonzero(self.database.selected.residues.resnames == lip))
            # sum_of_all_contacts += list(pivot_t.sum(axis=1) / n_frames)

        contact_metrics = pd.DataFrame({'Protein' : protein})
        contact_metrics['ResID'] = res_ids
        contact_metrics['ResName'] = res_names
        contact_metrics['Lipids'] = lipids
        contact_metrics['Radius'] = radius
        contact_metrics['Occupancy'] = occupancy
        contact_metrics['Lipid_Number'] = lipid_number
        # contact_metrics['Sum_of_all_Contacts'] = sum_of_all_contacts

        if save_file != '':
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
        protein = 'GIRK' # TODO: we'll need to update this into a list and iterate over it
        lipids = list(self.database_unique)
        sub_data = {k: {"category": k, "value": 0} for k in lipids} # TODO: we need to generate sub_data for each protein.
        js = {protein : {k: [] for k in lipids}}
        for residue, contact_counter in self.contacts_sum.items():
            for lipid, contact_sum in contact_counter.items():
                sub_data[lipid]['value'] += contact_sum
                metric = (contact_sum * self.dt) / self.totaltime
                js[protein][lipid].append([f'{resnames[residue]} {residue}', float("{:.2f}".format(metric))])

        sub_data = list(sub_data.values())
        norm_with = sum([x['value'] for x in sub_data])
        sub_data = [{'category': d['category'], 'value': "{:.2f}".format(d['value']/norm_with)} for d  in sub_data]

        # return js, {protein: sub_data}

        # Hardcoded
        proteins = ['GIRK']
        protein_counts = {'GIRK': 1}

        pie_data = []
        for protein in proteins:
            value = protein_counts[protein] / sum(protein_counts.values())

            protein_pdata = {
                "category": protein,
                "value": float("{:.2f}".format(value)),
                "subData": sub_data
            }
            pie_data.append(protein_pdata)

        # payload should include the entire data. The backend can process it then based on client requests
        payload = {
            "data": js,
            "proteins": [protein],
            "lipids": lipids,
            "pie_data": pie_data # TODO: include protein info
        }

        return payload


    def export_to_prolint(self, path='prolint_results.pkl'):
        """
        Temporal method to be able to use the analysis tools from **prolintpy**.

        Parameters
        ----------
        path : file
            Path to file to save the contacts information on Prolint's format. ('results_prolint.pkl')
            Prolint's results are stored as a dictionary. Contacts for each residue of each protein
            are stored using the ProLint.LPContacts class.
        """
        timestep = self.query.selected.universe.trajectory.dt

        PLASMA_LIPIDS = {}
        for lip in np.unique(self.database.selected.residues.resnames):
            PLASMA_LIPIDS[lip] = [lip]

        prolint_contacts = defaultdict(dict)
        n_residues_db = self.database.selected.residues.n_residues
        frames = self.database.selected.universe.trajectory.n_frames
        for protein in np.unique(self.query.selected.residues.macros):
            residues = self.query.selected.residues[self.query.selected.residues.macros == protein]
            per_residue_results = {}
            for idx in (pbar := tqdm(residues.resindices)):
                pbar.set_description('Exporting to Prolint in a per residue basis over protein residues')
                per_residue_results[idx+1] = LPContacts(self.contacts, self.counts, n_residues_db, frames, PLASMA_LIPIDS, self.database, timestep, residue=idx)

            prolint_contacts[protein][0] = per_residue_results

        with open(path, 'wb') as f:
            pickle.dump(prolint_contacts, f)

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


