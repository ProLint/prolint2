import contextlib
from collections import defaultdict
import joblib
from joblib import Parallel, delayed
from tqdm import tqdm
from MDAnalysis.analysis.base import AnalysisBase
from prolint2.computers.base import ContactComputerBase

class SerialContacts(AnalysisBase, ContactComputerBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class and the ContactComputerBase class.
    """
    def __init__(self, universe, query, database, cutoff, **kwargs):
        ContactComputerBase.__init__(self, universe, query, database, cutoff, **kwargs)
        AnalysisBase.__init__(self, self._universe.trajectory, **kwargs)

    def _single_frame(self):
        """
        Compute the contacts for a single frame.
        """
        pairs = self._compute_pairs()
        self._analyze_pairs_and_update(pairs, self._frame_index)

    def _analyze_pairs_and_update(self, pairs, frame_index):
        """
        Analyze the pairs of residues and lipids that are within the cutoff distance.
        """
        self._analyze_pairs(pairs, frame_index, self.contact_frames)

class ParallelContacts(ContactComputerBase):
    r"""

    Class to get the distance-based contacts starting from two AtomGroups
    using a *parallel* approach.

    It inherits from the ContactComputerBase class.
    """
    def __init__(self, universe, query, database, cutoff, n_blocks=4, **kwargs):
        self.n_blocks = n_blocks
        super().__init__(universe, query, database, cutoff, **kwargs)

    def analyze_block(self, blockslice):
        """
        Analyze a block of frames.
        """
        contact_frames = self._create_and_analyze_pairs(blockslice)
        return contact_frames

    def _create_and_analyze_pairs(self, blockslice):
        """
        Create and analyze the pairs of residues and lipids that are within the cutoff distance.
        """
        contact_frames = defaultdict(lambda: defaultdict(list))

        for ts in self._universe.trajectory[blockslice.start:blockslice.stop]:
            pairs = self._compute_pairs()
            self._analyze_pairs(pairs, ts.frame, contact_frames)
        
        return contact_frames

    def _get_block_slices(self):
        """
        Get the block slices.
        """
        n_frames = self._universe.trajectory.n_frames
        n_frames_per_block = n_frames // self.n_blocks
        blocks = [range(i * n_frames_per_block, (i + 1) * n_frames_per_block) for i in range(self.n_blocks - 1)]
        blocks.append(range((self.n_blocks - 1) * n_frames_per_block, n_frames))
        return blocks

    def run(self, **kwargs):
        """
        Run the analysis.
        """
        blocks = self._get_block_slices()
        # results = Parallel(n_jobs=self.n_blocks)(delayed(self.analyze_block)(bs) for bs in blocks)
        with tqdm_joblib(tqdm(total=self.n_blocks)) as _:
            results = Parallel(n_jobs=16)(delayed(self.analyze_block)(bs) for bs in blocks)

        # Merge results
        for contact_frames in results:
            for residue_id, lipid_data in contact_frames.items():
                for lipid_id, frames in lipid_data.items():
                    self.contact_frames[residue_id][lipid_id].extend(frames)
                    
        return self

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        """
        Class to patch joblib to report into tqdm progress bar given as argument.
        """
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
