"""Analysis module for ProLint.

Provides analysis classes for extracting insights from contact data:

- TimeSeriesAnalysis: Contact counts over trajectory time
- DatabaseContactsAnalysis: Per-molecule contact timelines
- KineticsAnalysis: Binding kinetics and residence times
- DensityMapAnalysis: 2D spatial density maps
- RadialDensityAnalysis: Radial density profiles
- DistanceAnalysis: Distance time series between residue pairs
- AtomDistancesAnalysis: Atom-atom distance matrices
- SharedContactsAnalysis: Pairwise residue correlations
- MetricsAnalysis: Per-residue contact metrics

All analyses are registered with AnalysisRegistry and can be accessed via
ComputedContacts.analyze() or directly instantiated.
"""

from prolint.analysis.base import BaseAnalysis, AnalysisRegistry, AnalysisResult
from prolint.analysis.density import DensityMapAnalysis, RadialDensityAnalysis
from prolint.analysis.shared_contacts import SharedContactsAnalysis
from prolint.analysis.timeseries import TimeSeriesAnalysis, DatabaseContactsAnalysis
from prolint.analysis.kinetics import KineticsAnalysis
from prolint.analysis.distances import DistanceAnalysis, AtomDistancesAnalysis
from prolint.analysis.metrics import MetricsAnalysis

# Register all built-in analyses
AnalysisRegistry.register("density_map", DensityMapAnalysis)
AnalysisRegistry.register("radial_density", RadialDensityAnalysis)
AnalysisRegistry.register("shared_contacts", SharedContactsAnalysis)
AnalysisRegistry.register("timeseries", TimeSeriesAnalysis)
AnalysisRegistry.register("database_contacts", DatabaseContactsAnalysis)
AnalysisRegistry.register("kinetics", KineticsAnalysis)
AnalysisRegistry.register("distances", DistanceAnalysis)
AnalysisRegistry.register("atom_distances", AtomDistancesAnalysis)
AnalysisRegistry.register("metrics", MetricsAnalysis)

__all__ = [
    "BaseAnalysis",
    "AnalysisRegistry",
    "AnalysisResult",
    "DensityMapAnalysis",
    "RadialDensityAnalysis",
    "SharedContactsAnalysis",
    "TimeSeriesAnalysis",
    "DatabaseContactsAnalysis",
    "KineticsAnalysis",
    "DistanceAnalysis",
    "AtomDistancesAnalysis",
    "MetricsAnalysis",
]
