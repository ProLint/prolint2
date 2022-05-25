"""
Ultr-Fast Contacts Calculation (UFCC)
=======================================

:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License

UFCC calculates de distance-based contacts between two references.
"""

import numpy as np
import warnings
import os
import scipy.stats
import scipy.sparse
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.analysis.base import AnalysisBase
from pmda.parallel import ParallelAnalysisBase
from multiprocessing import Pool

# import logging
# MDAnalysis.start_logging()

# logger = logging.getLogger("MDAnalysis.MDAKit.membrane_curvature")


class Contacts(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups.
    """

    def __init__(self, universe, query, database, 
                cutoff, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.cutoff = cutoff

        # to allow for non-sequential resindices
        self._sorted_protein_resindices = scipy.stats.rankdata(
            self.query.resindices,
            method="dense"
        ) - 1
        self._sorted_membrane_resindices = scipy.stats.rankdata(
            self.database.resindices,
            method="dense"
        ) - 1

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

        # # Apply wrapping coordinates
        # if not self.wrap:
        #     # Warning
        #     msg = (" `wrap == False` may result in inaccurate calculation "
        #            "of membrane curvature. Surfaces will be derived from "
        #            "a reduced number of atoms. \n "
        #            " Ignore this warning if your trajectory has "
        #            " rotational/translational fit rotations! ")
        #     warnings.warn(msg)
        #     logger.warn(msg)

    def _prepare(self):
        # Initialize empty np.array with results
        self.contacts = np.zeros(self.n_frames, dtype=object)

    def _single_frame(self):
        # Get the results and populate the results dictionary
        pairs = capped_distance(
            self.query.positions,
            self.database.positions,
            max_cutoff=self.cutoff,
            box=self.database.dimensions,
            return_distances=False
        )
        
        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        query_residx, database_residx = np.unique(np.array([[self._sorted_protein_resindices[pair[0]], self._sorted_membrane_resindices[pair[1]]] for pair in pairs]), axis=0).T

        # store neighbours for this frame
        data = np.ones_like(query_residx)
        self.contacts[self._frame_index] = scipy.sparse.csr_matrix((data, (query_residx, database_residx)), dtype=np.int8, shape=(self.query.n_residues, self.database.n_residues))
    
    # def _conclude(self):
    #     # OPTIONAL
    #     # Called once iteration on the trajectory is finished.
    #     # Apply normalisation and averaging to results here.
    #     self.result = np.asarray(self.result) / np.sum(self.result)
        

# class RMSD(ParallelAnalysisBase):
#     def __init__(self, mobile, ref, superposition=True):
#         universe = mobile.universe
#         super(RMSD, self).__init__(universe, (mobile, ))
#         self._ref_pos = ref.positions.copy()
#         self.superposition = superposition

#     def _prepare(self):
#         self.rmsd = None

#     def _conclude(self):
#         self.rmsd = np.vstack(self._results)

#     def _single_frame(self, ts, atomgroups):
#         return (ts.frame, ts.time,
#                 rms.rmsd(atomgroups[0].positions, self._ref_pos,
#                          superposition=self.superposition))
        
        
class ContactsPar(ParallelAnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups.
    """

    def __init__(self, universe, query, database, 
                cutoff, **kwargs):

        super().__init__(universe.universe.trajectory, (query, database))
        self.query = query
        self.query_n_residues = self.query.n_residues
        self.database = database
        self.database_n_residues = self.database.n_residues
        self.cutoff = cutoff

        # to allow for non-sequential resindices
        self._sorted_protein_resindices = scipy.stats.rankdata(
            self.query.resindices,
            method="dense"
        ) - 1
        self._sorted_membrane_resindices = scipy.stats.rankdata(
            self.database.resindices,
            method="dense"
        ) - 1

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

        # # Apply wrapping coordinates
        # if not self.wrap:
        #     # Warning
        #     msg = (" `wrap == False` may result in inaccurate calculation "
        #            "of membrane curvature. Surfaces will be derived from "
        #            "a reduced number of atoms. \n "
        #            " Ignore this warning if your trajectory has "
        #            " rotational/translational fit rotations! ")
        #     warnings.warn(msg)
        #     logger.warn(msg)

    def _prepare(self):
        # Initialize empty np.array with results
        self.contacts = None

    def _single_frame(self, ts, atomgroups):
        # Get the results and populate the results dictionary
        pairs = capped_distance(
            atomgroups[0].positions,
            atomgroups[1].positions,
            max_cutoff=self.cutoff,
            box=atomgroups[1].dimensions,
            return_distances=False
        )
        
        # Find unique pairs of residues interacting
        # Currently we have pairs of atoms
        query_residx, database_residx = np.unique(np.array([[self._sorted_protein_resindices[pair[0]], self._sorted_membrane_resindices[pair[1]]] for pair in pairs]), axis=0).T

        # store neighbours for this frame
        data = np.ones_like(query_residx)
        return (ts.frame, scipy.sparse.csr_matrix((data, (query_residx, database_residx)), dtype=np.int8, shape=(self.query_n_residues, self.database_n_residues)))
    
    def _conclude(self):
        self.contacts = np.array([l[1] for l in sorted(np.vstack(self._results), key=lambda tup: tup[0])])
        # self.contacts = self.contacts[:,1]
    #     # OPTIONAL
    #     # Called once iteration on the trajectory is finished.
    #     # Apply normalisation and averaging to results here.
    #     self.result = np.asarray(self.result) / np.sum(self.result)
               
        
        
        
        
        
        # gridsearch = FastNS(self.radius, self.database.positions, box=self.database.dimensions, pbc=True)
        # result = gridsearch.search(self.query.positions)
        # pairs = result.get_neighbours()
        # contacts_per_idx = []
        # flag = 0
        # for pair_idx in range(len(pairs)):
        #     if pair_idx == 0:
        #         contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
        #     else:
        #         if pairs[pair_idx][0] == pairs[pair_idx-1][0]:
        #             contacts_per_idx[flag].append(pairs[pair_idx][1])
        #         else:
        #             flag += 1
        #             contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
        # self.results.contacts[self._frame_index] = contacts_per_idx


class ContactsMP(object):
    """
    Base class for getting the contacts. 

    ...

    Atributes
    ---------
    u : MDA.Universe
        The whole system.
    query : MDA.AtomGroup 
        Query atoms.
    database : MDA.AtomGroup 
        database atoms.

    Methods
    -------
    get_contacts(n_jobs)
        Return the index pairs for the contacts on each frame.
    """

    def __init__(self, UFCC_universe, query, database, cutoff, n_jobs):
        """
        Parameters
        ---------
        u : MDA.Universe
            The whole system.
        query : MDA.AtomGroup 
            Query atoms.
        database : MDA.AtomGroup 
            Database atoms.
        """
        self.u = UFCC_universe
        self.query = query
        self.database = database
        self.n_jobs = n_jobs
        if self.n_jobs == -1:
            self.n_jobs = os.cpu_count()
        self.cutoff = cutoff
        
        # to allow for non-sequential resindices
        self._sorted_protein_resindices = scipy.stats.rankdata(
            self.query.resindices,
            method="dense"
        ) - 1
        self._sorted_membrane_resindices = scipy.stats.rankdata(
            self.database.resindices,
            method="dense"
        ) - 1

        self.contacts = None


    def per_block(self, blocksize):
        results = {}
        for ts in self.u.trajectory[blocksize.start:blocksize.stop]:
            # Get the results and populate the results dictionary
            pairs = capped_distance(
                self.query.positions,
                self.database.positions,
                max_cutoff=self.cutoff,
                box=self.database.dimensions,
                return_distances=False
            )
            
            # Find unique pairs of residues interacting
            # Currently we have pairs of atoms
            query_residx, database_residx = np.unique(np.array([[self._sorted_protein_resindices[pair[0]], self._sorted_membrane_resindices[pair[1]]] for pair in pairs]), axis=0).T

            # store neighbours for this frame
            data = np.ones_like(query_residx)
            results[ts.frame] = scipy.sparse.csr_matrix((data, (query_residx, database_residx)), dtype=np.int8, shape=(self.query.n_residues, self.database.n_residues))
        return results

    def make_blocks(self, n_blocks):
        n_frames = self.u.trajectory.n_frames  #len(prot)
        n_frames_per_block = n_frames // n_blocks
        blocks = [range(i * n_frames_per_block, (i + 1) * n_frames_per_block) for i in range(n_blocks - 1)]
        blocks.append(range((n_blocks - 1) * n_frames_per_block, n_frames))
        return blocks

    def run(self):
        #self.u.transfer_to_memory(verbose=True)
        blocks = self.make_blocks(self.n_jobs)
        with Pool(processes=self.n_jobs) as worker_pool:
            all_result = worker_pool.map(self.per_block, blocks)
            worker_pool.close()
        return all_result

    def get_contacts(self):
        """
        Parameters
        ---------
        get_contacts(n_jobs)
            Return the index pairs for the contacts on each frame.
        """
        dirty_contacts = self.run()
        clean_contacts = {}
        for block in dirty_contacts:
            clean_contacts = {**clean_contacts, **block}
        self.contacts = []
        for frame in sorted(clean_contacts.keys()):
            self.contacts.append(clean_contacts[frame])
        self.contacts = np.array(self.contacts)
