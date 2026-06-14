"""Distance analysis for residue pair interactions."""

from typing import Optional
import numpy as np
from MDAnalysis.lib.distances import distance_array

from prolint.analysis.base import BaseAnalysis, AnalysisResult


class AtomDistancesAnalysis(BaseAnalysis):
    """Compute atom-atom distance matrix at a specific frame.

    Provides detailed per-atom distance information between two residues
    at a single trajectory frame.

    See Also
    --------
    DistanceAnalysis : Distance time series over multiple frames
    """

    name = "atom_distances"
    """Analysis name for registry."""

    description = "Atom-atom distance matrix at a specific frame"
    """Human-readable description."""

    def run(
        self,
        query_residue: int,
        database_residue: int,
        frame_idx: int,
    ) -> AnalysisResult:
        """Compute atom-atom distance matrix between two residues.

        Parameters
        ----------
        query_residue : int
            Query residue ID.
        database_residue : int
            Database residue ID.
        frame_idx : int
            Frame index to analyze. Clamped to valid range.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - frame : int frame index
            - query_atoms : list of str atom names
            - database_atoms : list of str atom names
            - distance_matrix : 2D list of float distances
            - min_distance, max_distance : float extrema

        Raises
        ------
        ValueError
            If query or database residue not found in universe.
        """
        query_atoms = self.universe.query.select_atoms(f"resid {query_residue}")
        db_atoms = self.universe.database.select_atoms(f"resid {database_residue}")

        if len(query_atoms) == 0:
            raise ValueError(f"Query residue {query_residue} not found")
        if len(db_atoms) == 0:
            raise ValueError(f"Database residue {database_residue} not found")

        frame_idx = max(0, min(frame_idx, self.universe.trajectory.n_frames - 1))
        self.universe.trajectory[frame_idx]

        dist_matrix = distance_array(query_atoms.positions, db_atoms.positions)

        return AnalysisResult(
            data={
                "frame": frame_idx,
                "query_atoms": [a.name for a in query_atoms],
                "database_atoms": [a.name for a in db_atoms],
                "distance_matrix": dist_matrix.tolist(),
                "min_distance": float(np.min(dist_matrix)),
                "max_distance": float(np.max(dist_matrix)),
            },
            metadata={
                "query_residue": query_residue,
                "database_residue": database_residue,
                "n_query_atoms": len(query_atoms),
                "n_database_atoms": len(db_atoms),
            },
        )


class DistanceAnalysis(BaseAnalysis):
    """Analyze distances between query and database residues over time.

    Computes center-of-mass distances and optionally minimum atom-atom
    distances between a query-database residue pair across trajectory frames.

    See Also
    --------
    AtomDistancesAnalysis : Full distance matrix at single frame
    """

    name = "distances"
    """Analysis name for registry."""

    description = "Distance computations between residue pairs"
    """Human-readable description."""

    def run(
        self,
        query_residue: int,
        database_residue: int,
        frame_start: int = 0,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
        compute_min_distances: bool = True,
        compute_positions: bool = True,
    ) -> AnalysisResult:
        """Compute distance time series between a residue pair.

        Parameters
        ----------
        query_residue : int
            Query residue ID.
        database_residue : int
            Database residue ID.
        frame_start : int, default=0
            First frame to process.
        frame_end : int, optional
            Last frame (exclusive). Defaults to total frames.
        frame_step : int, default=1
            Step between frames.
        compute_min_distances : bool, default=True
            Whether to compute minimum atom-atom distances.
        compute_positions : bool, default=True
            Whether to compute 2D positions relative to query COM.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - frames : list of int frame indices
            - distances : list of float center-of-mass distances
            - contact_frames : list of int frames where contact occurred
            - min_distances : list of float (if compute_min_distances)
            - positions : dict with query/database 2D positions (if compute_positions)

        Raises
        ------
        ValueError
            If query or database residue not found in universe.
        """
        query_atoms = self.universe.query.select_atoms(f"resid {query_residue}")
        db_atoms = self.universe.database.select_atoms(f"resid {database_residue}")

        if len(query_atoms) == 0:
            raise ValueError(f"Query residue {query_residue} not found")
        if len(db_atoms) == 0:
            raise ValueError(f"Database residue {database_residue} not found")

        frames = self._get_frame_range(frame_start, frame_end, frame_step)

        # Get contact frames for this residue pair from contacts
        contact_frames_set = set()
        if query_residue in self.contacts.contact_frames:
            if database_residue in self.contacts.contact_frames[query_residue]:
                contact_frames_set = set(
                    self.contacts.contact_frames[query_residue][database_residue]
                )

        # Full query for computing relative positions
        full_query_atoms = self.universe.query.atoms

        com_distances = []
        min_distances = []
        contact_frames_list = []
        query_positions = []
        db_positions = []

        for frame_idx in frames:
            self.universe.trajectory[frame_idx]

            query_com = query_atoms.center_of_mass()
            db_com = db_atoms.center_of_mass()
            com_distances.append(float(np.linalg.norm(query_com - db_com)))

            if compute_positions:
                full_query_com = full_query_atoms.center_of_mass()[:2]
                query_positions.append(
                    {
                        "x": float(query_com[0] - full_query_com[0]),
                        "y": float(query_com[1] - full_query_com[1]),
                    }
                )
                db_positions.append(
                    {
                        "x": float(db_com[0] - full_query_com[0]),
                        "y": float(db_com[1] - full_query_com[1]),
                    }
                )

            if compute_min_distances:
                dist_matrix = distance_array(query_atoms.positions, db_atoms.positions)
                min_distances.append(float(np.min(dist_matrix)))

            if frame_idx in contact_frames_set:
                contact_frames_list.append(frame_idx)

        result = {
            "frames": frames,
            "distances": com_distances,
            "contact_frames": contact_frames_list,
        }

        if compute_min_distances:
            result["min_distances"] = min_distances

        if compute_positions:
            result["positions"] = {
                "query": query_positions,
                "database": db_positions,
            }

        return AnalysisResult(
            data=result,
            metadata={
                "query_residue": query_residue,
                "database_residue": database_residue,
                "frame_start": frame_start,
                "frame_end": (
                    frame_end if frame_end else self.universe.trajectory.n_frames
                ),
                "frame_step": frame_step,
                "n_frames": len(frames),
            },
        )
