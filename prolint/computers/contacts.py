"""Serial contact computation module.

This module provides the SerialContacts class for efficient
distance-based contact detection using grid search.
"""

import logging
import time
from collections import defaultdict

import numpy as np

from MDAnalysis.lib.nsgrid import FastNS

from prolint.computers.base import ContactComputerBase
from prolint.utils.utils import fast_unique_comparison

logger = logging.getLogger(__name__)


class SerialContacts(ContactComputerBase):
    """Distance-based contact detection using MDAnalysis FastNS.

    Computes contacts between query and database atom groups using
    a grid-based neighbor search algorithm for efficiency.

    Parameters
    ----------
    universe : Universe
        ProLint Universe instance.
    query : ExtendedAtomGroup
        Query atoms (e.g., protein).
    database : ExtendedAtomGroup
        Database atoms (e.g., lipids).
    cutoff : float
        Distance cutoff in Angstroms.
    **kwargs : dict
        Additional arguments passed to MDAnalysis AnalysisBase.

    Examples
    --------
    >>> from prolint import Universe
    >>> u = Universe("topology.gro", "trajectory.xtc")
    >>> contacts = u.compute_contacts(cutoff=7.0)

    See Also
    --------
    ContactComputerBase : Abstract base class
    ContactsProvider : Orchestrates contact computation
    """

    def __init__(self, universe, query, database, cutoff, **kwargs):
        super().__init__(universe.universe.trajectory, **kwargs)

        self.query = query
        self.database = database
        self.cutoff = cutoff
        self.contacts = None
        self.contact_frames = defaultdict(lambda: defaultdict(list))
        self._start_time = None
        self._last_log_frame = 0
        self._log_interval = 10  # Log every 10% progress

        self._validate_inputs()
        logger.debug("SerialContacts initialized with cutoff=%.1f Å", cutoff)

    def _validate_inputs(self):
        """Validate query and database selections.

        Raises
        ------
        ValueError
            If query or database is empty, or cutoff is non-positive.
        """
        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _compute_pairs(self):
        """Compute atom pairs within cutoff distance.

        Uses MDAnalysis FastNS grid search for efficient neighbor finding.

        Returns
        -------
        np.ndarray
            Array of shape (N, 2) containing query-database atom index pairs.
        """
        if self.database.dimensions is None:
            dim_x = np.max(self.database.universe.atoms.positions[:, 0]) - np.min(
                self.database.universe.atoms.positions[:, 0]
            )
            dim_y = np.max(self.database.universe.atoms.positions[:, 1]) - np.min(
                self.database.universe.atoms.positions[:, 1]
            )
            dim_z = np.max(self.database.universe.atoms.positions[:, 2]) - np.min(
                self.database.universe.atoms.positions[:, 2]
            )
            self.database.dimensions = np.array([dim_x, dim_y, dim_z, 90, 90, 90])

        gridsearch = FastNS(
            self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True
        )
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        return pairs

    def _prepare(self):
        """Called before iteration starts."""
        self._start_time = time.time()
        self._last_log_frame = 0
        logger.info(
            "Computing contacts: cutoff=%.1f Å, frames=%d",
            self.cutoff,
            self.n_frames,
        )

    def _single_frame(self):
        """Process a single trajectory frame.

        Computes contact pairs and stores frame indices for each
        residue-molecule contact.
        """
        pairs = self._compute_pairs()

        q_resid_indices = pairs[:, 0]
        db_resid_indices = pairs[:, 1]
        residue_ids = self.query.resids[q_resid_indices]
        database_ids = self.database.resids[db_resid_indices]
        database_names = self.database.resnames[db_resid_indices]

        residue_ids, database_ids, database_names = fast_unique_comparison(
            residue_ids, database_ids, database_names
        )

        existing_pairs = set()
        for unique_data in zip(residue_ids, database_ids, database_names):
            residue_id, database_id, _ = unique_data
            if (residue_id, database_id) not in existing_pairs:
                existing_pairs.add((residue_id, database_id))
                # NOTE: we need actual frame number (_ts.frame) instead of relative index (_frame_index) below
                # This ensures correct mapping when start != 0
                self.contact_frames[residue_id][database_id].append(self._ts.frame)

        # Log progress at intervals
        if self.n_frames > 0:
            progress = ((self._frame_index + 1) / self.n_frames) * 100
            if progress >= self._last_log_frame + self._log_interval:
                self._last_log_frame = (
                    int(progress // self._log_interval) * self._log_interval
                )
                logger.info(
                    "Frame progress: %d/%d (%.0f%%)",
                    self._frame_index + 1,
                    self.n_frames,
                    progress,
                )
        logger.debug("Frame %d: %d contact pairs", self._ts.frame, len(existing_pairs))

    def _conclude(self):
        """Called after iteration finishes."""
        elapsed = time.time() - self._start_time if self._start_time else 0
        total_contacts = sum(
            len(db_contacts)
            for residue_contacts in self.contact_frames.values()
            for db_contacts in residue_contacts.values()
        )
        logger.info(
            "Contact computation complete: %d residue-molecule pairs in %.2fs",
            total_contacts,
            elapsed,
        )
