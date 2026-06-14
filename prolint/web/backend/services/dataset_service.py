"""Dataset Service - handles loading and managing molecular dynamics datasets."""

from pathlib import Path
from typing import Optional
import logging
import tempfile
import uuid

from prolint import Universe
from prolint.core.replica_detection import detect_replicas
from prolint.web.backend.services.storage_service import storage

logger = logging.getLogger(__name__)

UPLOAD_DIR = Path(tempfile.gettempdir()) / "prolint_uploads"
UPLOAD_DIR.mkdir(exist_ok=True)


class DatasetService:
    """Service for loading and managing molecular dynamics datasets.

    Handles file uploads, Universe creation, and dataset storage
    for the web interface.
    """

    @staticmethod
    async def load_from_upload(
        topology_file,
        topology_filename: str,
        trajectory_file=None,
        trajectory_filename: Optional[str] = None,
        name: Optional[str] = None,
    ) -> str:
        """Load a dataset from uploaded files.

        Parameters
        ----------
        topology_file : file-like
            Topology file object (PDB, GRO, PSF, etc.).
        topology_filename : str
            Original filename for the topology.
        trajectory_file : file-like, optional
            Trajectory file object (XTC, TRR, DCD, etc.).
        trajectory_filename : str, optional
            Original filename for the trajectory.
        name : str, optional
            Custom name for the dataset.

        Returns
        -------
        str
            Unique dataset ID.

        Raises
        ------
        ValueError
            If files cannot be loaded as a valid Universe.
        """
        upload_id = str(uuid.uuid4())[:8]
        upload_subdir = UPLOAD_DIR / upload_id
        upload_subdir.mkdir(exist_ok=True)

        try:
            # Save topology
            topo_path = upload_subdir / topology_filename
            topology_file.seek(0)
            topo_path.write_bytes(topology_file.read())
            logger.info(f"Saved topology: {topo_path}")

            # Save trajectory if provided
            traj_path = None
            if trajectory_file and trajectory_filename:
                traj_path = upload_subdir / trajectory_filename
                trajectory_file.seek(0)
                traj_path.write_bytes(trajectory_file.read())
                logger.info(f"Saved trajectory: {traj_path}")

            # Load universe
            if traj_path:
                universe = Universe(str(topo_path), str(traj_path))
            else:
                universe = Universe(str(topo_path))

            logger.info(
                f"Loaded dataset: {universe.atoms.n_atoms} atoms, {universe.trajectory.n_frames} frames"
            )

            return storage.add_dataset(
                universe=universe,
                name=name or topo_path.stem,
                topology_file=str(topo_path),
                trajectory_file=str(traj_path) if traj_path else None,
            )

        except Exception as e:
            import shutil

            if upload_subdir.exists():
                shutil.rmtree(upload_subdir)
            raise ValueError(f"Could not load dataset: {e}")

    @staticmethod
    def get_dataset_info(dataset_id: str) -> Optional[dict]:
        """Get metadata for a specific dataset.

        Parameters
        ----------
        dataset_id : str
            Unique dataset identifier.

        Returns
        -------
        dict or None
            Dataset metadata if found, None otherwise.
        """
        return storage.get_dataset_metadata(dataset_id)

    @staticmethod
    def list_datasets() -> list[dict]:
        """List all loaded datasets.

        Returns
        -------
        list of dict
            Metadata for all datasets in storage.
        """
        return storage.list_datasets()

    @staticmethod
    def delete_dataset(dataset_id: str) -> bool:
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
        return storage.delete_dataset(dataset_id)

    @staticmethod
    def analyze_query_replicas(dataset_id: str, query_selection: str) -> dict:
        """Analyze replica structure in a query selection using core library.

        Detects if the query selection contains multiple replicas
        (e.g., protein replicates) and checks for repeated residue IDs.

        Parameters
        ----------
        dataset_id : str
            ID of the dataset to analyze.
        query_selection : str
            MDAnalysis selection string for query group.

        Returns
        -------
        dict
            Replica analysis results with keys:
            - success: bool
            - n_replicas: int
            - n_atoms: int
            - n_residues: int
            - detection_method: str or None
            - has_repeated_resids: bool or None
            - replica_info: list of replica details
            - message: str
        """
        universe = storage.get_dataset(dataset_id)
        if universe is None:
            raise ValueError(f"Dataset not found: {dataset_id}")

        # Select query atoms
        try:
            query_atoms = universe.select_atoms(query_selection)
        except Exception as e:
            raise ValueError(f"Invalid selection: {e}")

        if len(query_atoms) == 0:
            raise ValueError(f"Selection '{query_selection}' matched no atoms")

        n_atoms = len(query_atoms)
        n_residues = len(query_atoms.residues)

        # Use the core library for replica detection
        replica_result = detect_replicas(query_atoms)

        # Build message based on detection results
        if replica_result.n_replicas == 1:
            message = "Single replica detected in selection."
        else:
            method_display = (
                "bond connectivity" if replica_result.detection_method == "bond_connectivity"
                else "residue ID patterns"
            )
            if replica_result.has_repeated_resids:
                message = f"Detected {replica_result.n_replicas} replicas with repeated residue IDs (via {method_display})."
            else:
                message = f"Detected {replica_result.n_replicas} replicas with unique residue IDs (via {method_display})."

        return {
            "success": True,
            "n_atoms": n_atoms,
            "n_residues": n_residues,
            "n_replicas": replica_result.n_replicas,
            "detection_method": replica_result.detection_method,
            "has_repeated_resids": replica_result.has_repeated_resids,
            "replica_info": [
                {
                    "replica_id": info.replica_id,
                    "n_residues": info.n_residues,
                    "first_resid": info.first_resid,
                    "last_resid": info.last_resid,
                }
                for info in replica_result.replica_info
            ],
            "message": message,
        }
