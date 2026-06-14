"""Per-residue contact metrics analysis."""

from typing import Optional, List
from prolint.analysis.base import BaseAnalysis, AnalysisResult


class MetricsAnalysis(BaseAnalysis):
    """Per-residue contact metrics analysis.

    Computes per-residue contact metrics (occupancy, mean, max, sum) and
    returns values for all query residues in the universe.

    See Also
    --------
    ExactContacts.compute_metric : Underlying metric computation
    """

    name = "metrics"
    """Analysis name for registry."""

    description = "Per-residue contact metrics (occupancy, mean, max, sum)"
    """Human-readable description."""

    def run(self, **kwargs) -> AnalysisResult:
        """Compute per-residue contact metrics.

        Parameters
        ----------
        metric : str, default="occupancy"
            Metric to compute. One of "occupancy", "mean", "max", "sum".
        database_type : str, optional
            Filter by database residue name (e.g., "CHOL").
        query_residues : list of int, optional
            Specific query residues to include. If None, includes all.

        Returns
        -------
        AnalysisResult
            Result with data containing:

            - residues : list of dict with "resid" and "resname"
            - values : list of float metric values
            - metric : str metric name
            - database_type : str or None
        """
        # Extract parameters from kwargs
        metric: str = kwargs.get("metric", "occupancy")
        database_type: Optional[str] = kwargs.get("database_type")
        query_residues: Optional[List[int]] = kwargs.get("query_residues")

        # Compute raw metrics using existing method
        raw_metrics = self.contacts.compute_metric(metric, target_resname=database_type)

        # Get query residues from universe
        if query_residues is None:
            universe_residues = list(self.universe.query.residues)
        else:
            query_set = set(query_residues)
            universe_residues = [
                r for r in self.universe.query.residues if int(r.resid) in query_set
            ]

        # Build residue list and values
        residues = []
        values = []

        for res in universe_residues:
            resid = int(res.resid)
            resname = res.resname

            # Get value from raw metrics
            if resid in raw_metrics:
                if database_type and database_type in raw_metrics[resid]:
                    value = raw_metrics[resid][database_type].get("global", 0.0)
                elif database_type is None:
                    # Sum across all database types if no specific type requested
                    value = sum(
                        db_data.get("global", 0.0)
                        for db_data in raw_metrics[resid].values()
                    )
                else:
                    value = 0.0
            else:
                value = 0.0

            residues.append({"resid": resid, "resname": resname})
            values.append(float(value))

        return AnalysisResult(
            data={
                "residues": residues,
                "values": values,
                "metric": metric,
                "database_type": database_type,
            },
            metadata={
                "n_residues": len(residues),
                "min_value": min(values) if values else 0.0,
                "max_value": max(values) if values else 0.0,
            },
        )
