"""Density map analysis for spatial distribution of database molecules."""

import logging
from typing import Optional, List
import numpy as np

from prolint.analysis.base import BaseAnalysis, AnalysisResult

logger = logging.getLogger(__name__)


class RadialDensityAnalysis(BaseAnalysis):
    """Radial density profile analysis.

    Computes radially-averaged density from a 2D density map, useful for
    analyzing the radial distribution of database molecules around the query.

    See Also
    --------
    DensityMapAnalysis : Generates the 2D density map input
    """

    name = "radial_density"
    """Analysis name for registry."""

    description = "Radial density profile from 2D density map"
    """Human-readable description."""

    def run(
        self,
        density: List[List[float]],
        x_edges: List[float],
        y_edges: List[float],
        n_bins: int = 50,
    ) -> AnalysisResult:
        """Compute radial density profile from 2D density map.

        Parameters
        ----------
        density : list of list of float
            2D density array from DensityMapAnalysis.
        x_edges : list of float
            X bin edges from DensityMapAnalysis.
        y_edges : list of float
            Y bin edges from DensityMapAnalysis.
        n_bins : int, default=50
            Number of radial bins.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - r_centers : list of float radial bin centers
            - radial_density : list of float density values
            - r_max : float maximum radius
        """
        density_arr = np.asarray(density)

        x_centers = [(x_edges[i] + x_edges[i + 1]) / 2 for i in range(len(x_edges) - 1)]
        y_centers = [(y_edges[i] + y_edges[i + 1]) / 2 for i in range(len(y_edges) - 1)]

        X, Y = np.meshgrid(x_centers, y_centers)
        R = np.sqrt(X**2 + Y**2)

        r_max = float(R.max())
        r_bins = np.linspace(0, r_max, n_bins + 1)
        r_centers = ((r_bins[:-1] + r_bins[1:]) / 2).tolist()

        radial_density = []
        for i in range(len(r_bins) - 1):
            mask = (R >= r_bins[i]) & (R < r_bins[i + 1])
            if mask.any():
                radial_density.append(float(np.mean(density_arr.T[mask])))
            else:
                radial_density.append(0.0)

        return AnalysisResult(
            data={
                "r_centers": r_centers,
                "radial_density": radial_density,
                "r_max": r_max,
            },
            metadata={"n_bins": n_bins},
        )


class DensityMapAnalysis(BaseAnalysis):
    """Compute 2D spatial density maps of database molecules around query.

    Computes the 2D spatial distribution of database molecule positions
    relative to the query center of mass over trajectory frames.

    See Also
    --------
    RadialDensityAnalysis : Radially-averaged density from this output
    """

    name = "density_map"
    """Analysis name for registry."""

    description = "2D spatial density of database molecules around query"
    """Human-readable description."""

    def run(
        self,
        frame_start: int = 0,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
        bins: int = 50,
        database_types: Optional[List[str]] = None,
    ) -> AnalysisResult:
        """Compute 2D density map of database molecules.

        Parameters
        ----------
        frame_start : int, default=0
            First frame to process.
        frame_end : int, optional
            Last frame (exclusive). Defaults to total frames.
        frame_step : int, default=1
            Step between frames.
        bins : int, default=50
            Number of bins in each dimension.
        database_types : list of str, optional
            Database residue names to include (e.g., ["CHOL", "POPC"]).
            If None, includes all database atoms.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - density : 2D list of float density values
            - x_edges, y_edges : list of float bin edges
            - x_centers, y_centers : list of float bin centers
            - query_density : 2D list of float query atom density
        """
        frames = self._get_frame_range(frame_start, frame_end, frame_step)

        logger.info(
            "Computing density map: %d frames, %d bins",
            len(frames),
            bins,
        )

        query_atoms = self.universe.query.atoms
        if database_types is None:
            database_atoms = self.universe.database.atoms
        else:
            selection = "resname " + " ".join(database_types)
            database_atoms = self.universe.database.select_atoms(selection)

        density_result = self._compute_density(
            query_atoms, database_atoms, frames, bins
        )
        logger.debug("Density map computed: grid size %dx%d", bins, bins)

        return AnalysisResult(
            data=density_result,
            metadata={
                "frame_start": frame_start,
                "frame_end": (
                    frame_end if frame_end else self.universe.trajectory.n_frames
                ),
                "frame_step": frame_step,
                "bins": bins,
                "database_types": database_types,
                "n_frames": len(frames),
            },
        )

    def _compute_density(
        self, query_atoms, database_atoms, frames: List[int], bins: int
    ) -> dict:
        """Compute density histograms over trajectory frames.

        Applies periodic boundary conditions to wrap database positions
        to their nearest image relative to the query center of mass.
        """
        n_frames = len(frames)
        n_db_atoms = len(database_atoms)
        n_query_atoms = len(query_atoms)

        all_db_positions = np.empty((n_frames * n_db_atoms, 2), dtype=np.float32)
        all_query_positions = np.empty((n_frames * n_query_atoms, 2), dtype=np.float32)

        for i, frame_idx in enumerate(frames):
            self.universe.trajectory[frame_idx]
            query_com = query_atoms.center_of_mass()[:2]

            # Get box dimensions for PBC wrapping
            box = self.universe.dimensions
            box_xy = box[:2] if box is not None else None

            # Database positions relative to query COM
            db_rel = database_atoms.positions[:, :2] - query_com

            # Apply minimum image convention for periodic boundaries
            if box_xy is not None and box_xy[0] > 0 and box_xy[1] > 0:
                db_rel[:, 0] -= box_xy[0] * np.round(db_rel[:, 0] / box_xy[0])
                db_rel[:, 1] -= box_xy[1] * np.round(db_rel[:, 1] / box_xy[1])

            db_start = i * n_db_atoms
            all_db_positions[db_start : db_start + n_db_atoms] = db_rel

            # Query positions relative to query COM
            query_rel = query_atoms.positions[:, :2] - query_com

            # Apply minimum image convention for query atoms too
            if box_xy is not None and box_xy[0] > 0 and box_xy[1] > 0:
                query_rel[:, 0] -= box_xy[0] * np.round(query_rel[:, 0] / box_xy[0])
                query_rel[:, 1] -= box_xy[1] * np.round(query_rel[:, 1] / box_xy[1])

            query_start = i * n_query_atoms
            all_query_positions[query_start : query_start + n_query_atoms] = query_rel

        # Compute symmetric bin edges centered on query COM
        all_positions = np.vstack([all_db_positions, all_query_positions])
        max_range = (
            max(
                np.max(np.abs(all_positions[:, 0])), np.max(np.abs(all_positions[:, 1]))
            )
            * 1.05
        )

        bin_edges = np.linspace(-max_range, max_range, bins + 1)
        hist_range = [[-max_range, max_range], [-max_range, max_range]]

        density, x_edges, y_edges = np.histogram2d(
            all_db_positions[:, 0],
            all_db_positions[:, 1],
            bins=bin_edges,
            range=hist_range,
            density=True,
        )

        query_density, _, _ = np.histogram2d(
            all_query_positions[:, 0],
            all_query_positions[:, 1],
            bins=bin_edges,
            range=hist_range,
            density=True,
        )

        # Pre-compute bin centers for plotting
        x_centers = ((x_edges[:-1] + x_edges[1:]) / 2).tolist()
        y_centers = ((y_edges[:-1] + y_edges[1:]) / 2).tolist()

        return {
            "density": density.tolist(),
            "x_edges": x_edges.tolist(),
            "y_edges": y_edges.tolist(),
            "x_centers": x_centers,
            "y_centers": y_centers,
            "query_density": query_density.tolist(),
        }
