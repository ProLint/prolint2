"""Shared contacts analysis for pairwise residue correlations."""

from typing import Optional, List, Dict, Set

from prolint.analysis.base import BaseAnalysis, AnalysisResult


class SharedContactsAnalysis(BaseAnalysis):
    """Analyze pairwise correlations between query residues.

    Identifies query residue pairs that interact with the same database
    molecule simultaneously, revealing potential cooperative binding sites.

    See Also
    --------
    TimeSeriesAnalysis : Contact dynamics over time
    """

    name = "shared_contacts"
    """Analysis name for registry."""

    description = (
        "Pairwise correlations between query residues via shared database contacts"
    )
    """Human-readable description."""

    def run(
        self,
        database_type: Optional[str] = None,
        normalize: bool = False,
    ) -> AnalysisResult:
        """Compute shared contacts correlation matrix.

        Parameters
        ----------
        database_type : str, optional
            Filter by database residue name (e.g., "CHOL").
            If None, includes all database molecules.
        normalize : bool, default=False
            Whether to normalize the matrix to 0-1 range.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - labels : list of int query residue IDs
            - matrix : 2D list of int/float shared contact counts
            - residue_to_index : dict mapping residue ID to matrix index
        """
        # Filter contacts by database type and build inverted index
        filtered_contacts = self._filter_by_database_type(database_type)
        database_to_query_frames = self._build_inverted_index(filtered_contacts)

        # Get all query residues from the inverted index
        all_query_residues: Set[int] = set()
        for query_frames_dict in database_to_query_frames.values():
            all_query_residues.update(query_frames_dict.keys())

        query_residue_list = sorted(int(r) for r in all_query_residues)
        n = len(query_residue_list)

        if n < 2:
            return AnalysisResult(
                data={
                    "labels": query_residue_list,
                    "matrix": [],
                    "residue_to_index": {},
                },
                metadata={
                    "database_type": database_type,
                    "n_residues": n,
                    "normalize": normalize,
                },
            )

        residue_to_index = {res: idx for idx, res in enumerate(query_residue_list)}
        matrix = [[0] * n for _ in range(n)]

        # Compute pairwise shared contacts
        for query_frames_dict in database_to_query_frames.values():
            query_resids = list(query_frames_dict.keys())
            for i in range(len(query_resids)):
                for j in range(i + 1, len(query_resids)):
                    res1, res2 = query_resids[i], query_resids[j]
                    shared_frames = len(
                        query_frames_dict[res1] & query_frames_dict[res2]
                    )
                    if shared_frames > 0:
                        idx1, idx2 = (
                            residue_to_index[int(res1)],
                            residue_to_index[int(res2)],
                        )
                        matrix[idx1][idx2] += shared_frames
                        matrix[idx2][idx1] += shared_frames

        if normalize and matrix:
            max_val = max(max(row) for row in matrix)
            if max_val > 0:
                matrix = [[v / max_val for v in row] for row in matrix]

        return AnalysisResult(
            data={
                "labels": query_residue_list,
                "matrix": matrix,
                "residue_to_index": residue_to_index,
            },
            metadata={
                "database_type": database_type,
                "n_residues": n,
                "normalize": normalize,
                "n_database_molecules": len(database_to_query_frames),
            },
        )

    def _build_inverted_index(
        self, filtered_contacts: Dict[int, Dict[int, List[int]]]
    ) -> Dict[int, Dict[int, Set[int]]]:
        """Build inverted index from contact frames.

        Parameters
        ----------
        filtered_contacts : dict
            Contact frames dict mapping query_resid -> database_id -> frame list.

        Returns
        -------
        dict
            Inverted mapping: database_id -> {query_resid: set of frames}.
        """
        database_to_query_frames: Dict[int, Dict[int, Set[int]]] = {}

        for query_resid, db_dict in filtered_contacts.items():
            for db_id, frames in db_dict.items():
                if db_id not in database_to_query_frames:
                    database_to_query_frames[db_id] = {}
                database_to_query_frames[db_id][query_resid] = set(frames)

        return database_to_query_frames
