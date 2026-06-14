"""Interaction computation service.

This module provides the service layer for computing biomolecular
contacts between query and database selections.
"""

from typing import Optional
import time
import logging

from prolint.web.backend.services.storage_service import storage, InteractionResult

logger = logging.getLogger(__name__)


class InteractionService:
    """Service for computing biomolecular interactions.

    Orchestrates contact computation using stored datasets and
    manages result storage for dashboard visualization.
    """

    @staticmethod
    def compute_all_granularities(
        dataset_id: str,
        groupA_selection: str,
        groupB_selection: str,
        cutoff: float = 7.0,
        start: int = 0,
        stop: Optional[int] = None,
        step: int = 1,
        units: str = "ns",
        normalize_by: str = "counts",
        selected_replica: Optional[str] = None,
        replica_info: Optional[list] = None,
    ) -> tuple[str, float]:
        """Compute contacts between two atom selections.

        Parameters
        ----------
        dataset_id : str
            ID of the dataset to analyze.
        groupA_selection : str
            MDAnalysis selection string for query group.
        groupB_selection : str
            MDAnalysis selection string for database group.
        cutoff : float, default=7.0
            Contact distance cutoff in Angstroms.
        start : int, default=0
            Starting frame index.
        stop : int, optional
            Ending frame index (None for last frame).
        step : int, default=1
            Frame step interval.
        units : str, default="ns"
            Time units for normalization.
        normalize_by : str, default="counts"
            Normalization method.
        selected_replica : str, optional
            Replica ID to analyze (e.g., 'A', 'B'). Required when multiple
            replicas with repeated residue IDs are detected.
        replica_info : list, optional
            Information about detected replicas for ViewPage display.

        Returns
        -------
        tuple of (str, float)
            Result ID and computation time in seconds.

        Raises
        ------
        ValueError
            If dataset not found or if repeated residue IDs detected
            without specifying selected_replica.
        """
        universe = storage.get_dataset(dataset_id)
        if universe is None:
            raise ValueError(f"Dataset not found: {dataset_id}")

        universe.units = units
        universe.normalize_by = normalize_by

        # Update query and database selections with validation
        try:
            universe.query = universe.select_atoms(groupA_selection)
        except Exception as e:
            raise ValueError(f"Invalid query selection '{groupA_selection}': {e}")

        if universe.query.n_atoms == 0:
            raise ValueError(f"Query selection '{groupA_selection}' matched no atoms")

        try:
            universe.database = universe.select_atoms(groupB_selection)
        except Exception as e:
            raise ValueError(f"Invalid database selection '{groupB_selection}': {e}")

        if universe.database.n_atoms == 0:
            raise ValueError(f"Database selection '{groupB_selection}' matched no atoms")

        logger.info(
            f"Computing: {groupA_selection} ({universe.query.n_atoms} atoms) <-> "
            f"{groupB_selection} ({universe.database.n_atoms} atoms), cutoff={cutoff}"
        )

        start_time = time.time()
        contacts = universe.compute_contacts(
            cutoff=cutoff,
            start=start,
            stop=stop,
            step=step,
            replica=selected_replica,
        )
        computation_time = time.time() - start_time

        n_frames = universe.trajectory.n_frames
        frame_end = min(stop, n_frames) if stop is not None else n_frames

        result = InteractionResult(
            universe=universe,
            contacts=contacts,
            computation_time=computation_time,
            frame_start=start,
            frame_end=frame_end,
            frame_step=step,
            units=units,
            normalize_by=normalize_by,
            norm_factor=universe.params["norm_factor"],
            selected_replica=selected_replica,
            replica_info=replica_info,
        )

        result_id = storage.add_result(result)
        logger.info(
            f"Computed in {computation_time:.2f}s "
            f"(selected replica: {selected_replica or 'all'})"
        )

        return result_id, computation_time
