"""
In-Memory Storage Service

Manages datasets, universes, and interaction results in memory.
"""

from typing import Dict, Optional, Any
import uuid
import logging
import numpy as np

from prolint import Universe


def _convert_numpy_types(obj: Any) -> Any:
    """Recursively convert numpy types to native Python types.

    Parameters
    ----------
    obj : any
        Object potentially containing numpy types.

    Returns
    -------
    any
        Object with numpy types converted to Python natives
        for JSON serialization.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, dict):
        return {
            _convert_numpy_types(k): _convert_numpy_types(v) for k, v in obj.items()
        }
    elif isinstance(obj, (list, tuple)):
        return [_convert_numpy_types(item) for item in obj]
    return obj


logger = logging.getLogger(__name__)


class InteractionResult:
    """Container for computed contacts and visualization metadata.

    Stores the Universe, ComputedContacts, and computation parameters
    needed for dashboard visualization endpoints.

    Attributes
    ----------
    universe : Universe
        The molecular dynamics Universe.
    contacts : ComputedContacts
        Computed contact data.
    computation_time : float
        Time taken for computation in seconds.
    frame_start : int
        Starting frame index.
    frame_end : int
        Ending frame index.
    frame_step : int
        Frame step interval.
    units : str
        Time units for display.
    normalize_by : str
        Normalization method used.
    norm_factor : float
        Normalization factor applied.
    selected_replica : str or None
        Selected replica for analysis, or None if all replicas.
    replica_info : list
        Information about detected replicas.
    """

    def __init__(
        self,
        universe: Universe,
        contacts,
        computation_time: float,
        frame_start: int = 0,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
        units: str = "ns",
        normalize_by: str = "actual_time",
        norm_factor: float = 1.0,
        selected_replica: Optional[str] = None,
        replica_info: Optional[list] = None,
    ):
        self.universe = universe
        self.contacts = contacts
        self.computation_time = computation_time
        self.frame_start = frame_start
        n_frames = universe.trajectory.n_frames
        self.frame_end = min(frame_end, n_frames) if frame_end is not None else n_frames
        self.frame_step = frame_step
        self.units = units
        self.normalize_by = normalize_by
        self.norm_factor = norm_factor
        self.selected_replica = selected_replica
        self.replica_info = replica_info or []

    def get_interaction_data(self) -> Dict[str, Any]:
        """Get interaction data for visualization.

        Returns
        -------
        dict
            Composition, universe info, frame range, and parameters.
        """
        return _convert_numpy_types(
            {
                "composition": {
                    "resname_counts": self.universe.database.resname_counts,
                },
                "universe": {
                    "n_frames": self.universe.trajectory.n_frames,
                },
                "frame_range": {
                    "start": self.frame_start,
                    "end": self.frame_end,
                    "step": self.frame_step,
                },
                "params": {
                    "units": self.units,
                    "normalize_by": self.normalize_by,
                    "norm_factor": self.norm_factor,
                    "selected_replica": self.selected_replica,
                    "replica_info": self.replica_info,
                },
            }
        )


class StorageService:
    """In-memory storage for datasets and computation results.

    Manages Universe instances and InteractionResult objects
    for the web dashboard.

    Attributes
    ----------
    _datasets : dict
        Mapping of dataset IDs to Universe instances.
    _results : dict
        Mapping of result IDs to InteractionResult instances.
    _dataset_metadata : dict
        Mapping of dataset IDs to metadata dicts.
    """

    def __init__(self):
        self._datasets: Dict[str, Universe] = {}
        self._results: Dict[str, InteractionResult] = {}
        self._dataset_metadata: Dict[str, dict] = {}

    # Dataset Management

    def add_dataset(
        self,
        universe: Universe,
        name: str,
        topology_file: Optional[str] = None,
        trajectory_file: Optional[str] = None,
    ) -> str:
        """Add a dataset to storage.

        Parameters
        ----------
        universe : Universe
            Loaded molecular dynamics Universe.
        name : str
            Display name for the dataset.
        topology_file : str, optional
            Path to topology file.
        trajectory_file : str, optional
            Path to trajectory file.

        Returns
        -------
        str
            Unique dataset ID.
        """
        dataset_id = str(uuid.uuid4())
        self._datasets[dataset_id] = universe
        self._dataset_metadata[dataset_id] = {
            "id": dataset_id,
            "name": name,
            "topology_file": topology_file,
            "trajectory_file": trajectory_file,
            "n_frames": universe.trajectory.n_frames,
            "n_atoms": universe.atoms.n_atoms,
            "n_residues": universe.atoms.n_residues,
            "status": "ready",
        }
        logger.info(
            f"Added dataset {dataset_id}: {name} ({universe.trajectory.n_frames} frames)"
        )
        return dataset_id

    def get_dataset(self, dataset_id: str) -> Optional[Universe]:
        """Get a dataset by ID.

        Parameters
        ----------
        dataset_id : str
            Unique dataset identifier.

        Returns
        -------
        Universe or None
            The Universe if found, None otherwise.
        """
        return self._datasets.get(dataset_id)

    def get_dataset_metadata(self, dataset_id: str) -> Optional[dict]:
        """Get metadata for a dataset.

        Parameters
        ----------
        dataset_id : str
            Unique dataset identifier.

        Returns
        -------
        dict or None
            Dataset metadata if found, None otherwise.
        """
        return self._dataset_metadata.get(dataset_id)

    def list_datasets(self) -> list[dict]:
        """List all datasets.

        Returns
        -------
        list of dict
            Metadata for all stored datasets.
        """
        return list(self._dataset_metadata.values())

    def delete_dataset(self, dataset_id: str) -> bool:
        """Delete a dataset from storage.

        Parameters
        ----------
        dataset_id : str
            Unique dataset identifier.

        Returns
        -------
        bool
            True if deleted, False if not found.
        """
        if dataset_id in self._datasets:
            del self._datasets[dataset_id]
            del self._dataset_metadata[dataset_id]
            logger.info(f"Deleted dataset {dataset_id}")
            return True
        return False

    # Result Management

    def add_result(self, result: InteractionResult) -> str:
        """Add a computation result to storage.

        Parameters
        ----------
        result : InteractionResult
            Computed interaction result.

        Returns
        -------
        str
            Unique result ID.
        """
        result_id = str(uuid.uuid4())
        self._results[result_id] = result
        logger.info(f"Added result {result_id}")
        return result_id

    def get_result(self, result_id: str) -> Optional[InteractionResult]:
        """Get a computation result by ID.

        Parameters
        ----------
        result_id : str
            Unique result identifier.

        Returns
        -------
        InteractionResult or None
            The result if found, None otherwise.
        """
        return self._results.get(result_id)


# Global storage instance
storage = StorageService()
