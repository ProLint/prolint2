"""Time series analysis for contact dynamics over trajectory."""

import logging
from typing import Optional, List

from prolint.analysis.base import BaseAnalysis, AnalysisResult

logger = logging.getLogger(__name__)


class DatabaseContactsAnalysis(BaseAnalysis):
    """Compute per-database-molecule contact timeline for a query residue.

    Creates a binary contact matrix showing which database molecules are
    in contact with a specific query residue at each frame.

    See Also
    --------
    TimeSeriesAnalysis : Aggregated contact counts over time
    """

    name = "database_contacts"
    """Analysis name for registry."""

    description = "Per-database-molecule contact timeline"
    """Human-readable description."""

    def run(
        self,
        query_residue: int,
        database_type: Optional[str] = None,
        frame_start: int = 0,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
        top_n: Optional[int] = None,
    ) -> AnalysisResult:
        """Compute binary contact matrix for a query residue.

        Parameters
        ----------
        query_residue : int
            Query residue ID to analyze.
        database_type : str, optional
            Filter by database residue name (e.g., "CHOL").
        frame_start : int, default=0
            First frame to process.
        frame_end : int, optional
            Last frame (exclusive). Defaults to total frames.
        frame_step : int, default=1
            Step between frames.
        top_n : int, optional
            If specified, return only the top N database IDs with the most
            contacts. Useful for reducing response size when many database
            molecules have minimal contacts.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - database_ids : list of int sorted database molecule IDs
            - frames : list of int frame indices
            - contact_matrix : dict mapping database_id to list of int (0 or 1)
        """
        frames = self._get_frame_range(frame_start, frame_end, frame_step)

        filtered_contacts = self._filter_by_database_type(database_type)
        db_dict = filtered_contacts.get(query_residue, {})

        database_ids = []
        contact_matrix = {}

        for db_id, contact_frames_list in db_dict.items():
            database_ids.append(int(db_id))
            contact_frame_set = set(contact_frames_list)
            contact_matrix[int(db_id)] = [
                1 if frame in contact_frame_set else 0 for frame in frames
            ]

        # Track total count before any filtering
        total_database_ids = len(database_ids)

        # Filter to top N database IDs by total contacts if specified
        if top_n is not None and len(database_ids) > top_n:
            # Calculate total contacts for each database_id
            contact_totals = {
                db_id: sum(contact_matrix[db_id]) for db_id in database_ids
            }
            # Sort by total contacts descending, then by ID for stability
            sorted_ids = sorted(
                database_ids, key=lambda x: (-contact_totals[x], x)
            )
            # Keep only top N
            top_ids = set(sorted_ids[:top_n])
            database_ids = [db_id for db_id in database_ids if db_id in top_ids]
            contact_matrix = {
                db_id: contact_matrix[db_id] for db_id in database_ids
            }

        return AnalysisResult(
            data={
                "database_ids": sorted(database_ids),
                "frames": frames,
                "contact_matrix": contact_matrix,
                "total_database_ids": total_database_ids,
            },
            metadata={
                "query_residue": query_residue,
                "database_type": database_type,
                "frame_start": frame_start,
                "frame_end": (
                    frame_end if frame_end else self.universe.trajectory.n_frames
                ),
                "frame_step": frame_step,
                "n_frames": len(frames),
                "top_n": top_n,
            },
        )


class TimeSeriesAnalysis(BaseAnalysis):
    """Analyze contact dynamics over trajectory time.

    Computes per-frame contact counts for query residues, useful for
    understanding how contact behavior varies throughout the simulation.

    See Also
    --------
    DatabaseContactsAnalysis : Per-molecule binary contact timeline
    KineticsAnalysis : Binding kinetics and residence times
    """

    name = "timeseries"
    """Analysis name for registry."""

    description = "Contact counts over trajectory time"
    """Human-readable description."""

    def run(
        self,
        database_type: Optional[str] = None,
        query_residues: Optional[List[int]] = None,
        frame_start: int = 0,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
    ) -> AnalysisResult:
        """Compute contact count time series.

        Parameters
        ----------
        database_type : str, optional
            Filter by database residue name (e.g., "CHOL").
        query_residues : list of int, optional
            Specific query residues to include. If None, includes all
            residues with contacts.
        frame_start : int, default=0
            First frame to process.
        frame_end : int, optional
            Last frame (exclusive). Defaults to total frames.
        frame_step : int, default=1
            Step between frames.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - query_residues : list of int query residue IDs
            - frames : list of int frame indices
            - contact_counts : dict mapping residue_id to list of counts
        """
        frames = self._get_frame_range(frame_start, frame_end, frame_step)
        n_frames = len(frames)
        actual_frame_end = frame_end if frame_end else self.universe.trajectory.n_frames

        filtered_contacts = self._filter_by_database_type(database_type)

        if query_residues is None:
            query_residues = sorted(filtered_contacts.keys())
        else:
            query_residues = [int(r) for r in query_residues]

        logger.info(
            "Computing time series: %d frames, %d residues",
            n_frames,
            len(query_residues),
        )

        contact_counts = {}
        for query_resid in query_residues:
            frame_counts = [0] * n_frames
            db_dict = filtered_contacts.get(query_resid, {})

            for contact_frames_list in db_dict.values():
                for frame in contact_frames_list:
                    if frame_start <= frame < actual_frame_end:
                        frame_idx = (frame - frame_start) // frame_step
                        if 0 <= frame_idx < n_frames:
                            frame_counts[frame_idx] += 1

            contact_counts[query_resid] = frame_counts

        return AnalysisResult(
            data={
                "query_residues": query_residues,
                "frames": frames,
                "contact_counts": contact_counts,
            },
            metadata={
                "database_type": database_type,
                "frame_start": frame_start,
                "frame_end": actual_frame_end,
                "frame_step": frame_step,
                "n_frames": n_frames,
            },
        )
