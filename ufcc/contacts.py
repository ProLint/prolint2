# Daniel P. Ramirez & Besian I. Sejdiu
# Prolint: A tool to analyze and visualize lipid-protein interactions.
#


import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS
from multiprocessing import Pool


class Contacts(object):
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

    def __init__(self, UFCC_universe, query, database):
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

    def per_block(self, blocksize):
        results = {}
        for ts in self.u.trajectory[blocksize.start:blocksize.stop]:
            gridsearch = FastNS(7, self.database.positions, box=self.database.dimensions, pbc=True)
            result = gridsearch.search(self.query.positions)
            pairs = result.get_pairs()
            contacts_per_idx = []
            flag = 0
            for pair_idx in range(len(pairs)):
                if pair_idx == 0:
                    contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
                else:
                    if pairs[pair_idx][0] == pairs[pair_idx-1][0]:
                        contacts_per_idx[flag].append(pairs[pair_idx][1])
                    else:
                        flag += 1
                        contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
            results[ts.frame] = contacts_per_idx
        return results

    def make_blocks(self, n_blocks):
        n_frames = self.u.trajectory.n_frames  #len(prot)
        n_frames_per_block = n_frames // n_blocks
        blocks = [range(i * n_frames_per_block, (i + 1) * n_frames_per_block) for i in range(n_blocks - 1)]
        blocks.append(range((n_blocks - 1) * n_frames_per_block, n_frames))
        return blocks

    def run(self, n_jobs):
        #self.u.transfer_to_memory(verbose=True)
        blocks = self.make_blocks(n_jobs)
        with Pool(processes=n_jobs) as worker_pool:
            all_result = worker_pool.map(self.per_block, blocks)
            worker_pool.close()
        return all_result

    def get_contacts(self, n_jobs):
        """
        Parameters
        ---------
        get_contacts(n_jobs)
            Return the index pairs for the contacts on each frame.
        """
        dirty_contacts = self.run(n_jobs)
        clean_contacts = {}
        for block in dirty_contacts:
            clean_contacts = {**clean_contacts, **block}
        return clean_contacts
        # TODO: create tool to process results
