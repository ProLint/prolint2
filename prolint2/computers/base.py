from collections import defaultdict
from MDAnalysis.lib.nsgrid import FastNS
from prolint2.utils.utils import fast_unique_comparison

class ContactComputerBase:
    """
    Base class for the contact computers.
    """
    def __init__(self, universe, query, database, cutoff, **kwargs):
        self._universe = universe.universe
        self.query = query
        self.database = database
        self.cutoff = cutoff
        self.q_resids = self.query.resids
        self.db_resids = self.database.resids
        self.db_resnames = self.database.resnames
        self.contacts = None
        self.contact_frames = defaultdict(lambda: defaultdict(list))
        self._validate_inputs()

    def _validate_inputs(self):
        """
        Validate inputs.
        """
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _get_residue_lipid_info(self, pair):
        """
        Get the residue and lipid information for a given pair.
        """
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

    def _analyze_pairs(self, pairs, frame_index, contact_frames):
        """
        Analyze the pairs of residues and lipids that are within the cutoff distance.
        """
        q_resid_indices = pairs[:, 0]
        db_resid_indices = pairs[:, 1]
        residue_ids = self.q_resids[q_resid_indices]
        lipid_ids = self.db_resids[db_resid_indices]
        lipid_names = self.db_resnames[db_resid_indices]

        residue_ids, lipid_ids, lipid_names = fast_unique_comparison(residue_ids, lipid_ids, lipid_names)

        existing_pairs = set()
        for unique_data in zip(residue_ids, lipid_ids, lipid_names):
            residue_id, lipid_id, _ = unique_data
            if (residue_id, lipid_id) not in existing_pairs:
                existing_pairs.add((residue_id, lipid_id))
                contact_frames[residue_id][lipid_id].append(frame_index)
