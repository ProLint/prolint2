r""":mod:`prolint2.metrics.base`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from abc import ABC, abstractmethod
from collections import defaultdict

from typing import Type, List, Union, Callable

from scipy.optimize import curve_fit
from prolint2.metrics.formatters import OutputFormat, DefaultOutputFormat


MetricRegistry = Type["registries.MetricRegistry"]


class BaseMetric(ABC):
    """Base class for all metrics classes that act on single frame contact Iterables."""

    name: str = None

    def __init__(self):
        pass

    @abstractmethod
    def compute_metric(self, contact_array):
        pass

    @classmethod
    def _register(cls, registry: MetricRegistry):
        registry.register(cls.name, cls)


class Metric(ABC):
    """Base class for metric calculation."""

    def __init__(
        self,
        contacts,
        metrics,
        output_format: OutputFormat = DefaultOutputFormat(),
        lipid_type=None,
        clear=True,
    ):
        self.contact_input = dict(sorted(contacts.contacts.items()))

        if not isinstance(metrics, list):
            metrics = [metrics]
        self.metrics = metrics
        if clear:
            output_format.clear()
        self.output_format = output_format
        self.lipid_type = lipid_type

    def compute(self, dt=1, totaltime=1):
        """Compute the metric for the given contacts."""
        multiplier = dt / totaltime
        for residue_id, lipid_dict in self.contact_input.items():
            for lipid_name, contact_array in lipid_dict.items():
                if self.lipid_type is not None and self.lipid_type != lipid_name:
                    continue
                # contact_array = list(lipid_contacts.values())

                if contact_array:
                    for metric in self.metrics:
                        # if max(contact_array) > 1:
                        #     print ('contact_array', residue_id, lipid_name, max(contact_array))
                        value = metric.compute_metric(contact_array) * multiplier
                        # print ('value', residue_id, lipid_name, value, multiplier)
                        self.output_format.store_result(
                            residue_id, lipid_name, metric.__class__.__name__, value
                        )
                else:
                    for metric in self.metrics:
                        self.output_format.store_result(
                            residue_id, lipid_name, metric.__class__.__name__, 0
                        )

        return self.output_format.get_result()


class BaseContactStore:
    """Base class for storing contact."""

    def __init__(self, ts, contact_frames, norm_factor: float = 1.0):
        self.norm_factor = float(norm_factor)
        self.contact_frames = contact_frames

        self._resids = ts.database.residues.resids
        self._resnames = ts.database.residues.resnames
        self._database_unique_resnames = ts.database.unique_resnames
        self._contacts = defaultdict(lambda: defaultdict(dict))
        self._nframes = ts.trajectory.n_frames

    def run(self, lipid_resnames: Union[str, List] = None):
        """Run the contact calculation for the given lipid resnames. If no resnames are given, all resnames are used."""
        raise NotImplementedError("Subclasses should implement this method.")

    def compute(self, metric: str, target_lipid_name=None):
        """Compute a pre-defined metric for all lipids or a specific lipid.

        Parameters
        ----------
        metric : str
            The metric to compute. Must be one of 'max', 'sum', 'mean'.
        target_lipid_name : str, optional
            The name of the lipid to compute the metric for. If None, the metric will be computed for all lipids.

        Returns
        -------
        Dict[str, Dict[str, Dict[int, float]]]
            A dictionary of computed metrics for all lipids.

        Examples
        --------
        >>> cd = AproxContacts(...)
        >>> cd.run()
        >>> cd.compute('max')
        >>> cd.compute('sum', 'DOPC')
        >>> cd.compute('median') # raises ValueError. Use `apply_function` instead.
        """

        if metric in ["max", "sum", "mean"]:
            return self.compute_metric(metric, target_lipid_name)
        else:
            raise ValueError(
                "Invalid metric specified. Use 'max', 'sum', 'mean'. For more complex metrics, use `apply_function`."
            )

    def compute_metric(self, metric: str, target_lipid_name=None):
        """Compute the given metric for the given lipid name."""
        raise NotImplementedError("Subclasses should implement this method.")

    def apply_function(self, func: Callable, target_lipid_name=None):
        """Apply the given function to the contacts for the given lipid name."""
        raise NotImplementedError("Subclasses should implement this method.")

    def pooled_results(self):
        """Get the computed contacts all pooled together."""
        raise NotImplementedError("Subclasses should implement this method.")

    @property
    def results(self):
        """Get the computed contacts per lipid id."""
        if self._contacts is None:
            raise ValueError("No contacts have been computed yet. Call run() first.")
        return self._contacts

    @property
    def contacts(self):
        """Get the computed contacts all pooled together."""
        if self._contacts is None:
            raise ValueError("No contacts have been computed yet. Call run() first.")
        return self.pooled_results()
    
    @property
    def occupancies(self):
        """Get the computed occupancies all pooled together."""
        if self.contact_frames is None:
            raise ValueError("No contacts have been computed yet. Call run() first.")
        occupancies_dict = {}
        for res_id in self.contact_frames.keys():
            occupancies_dict[res_id] = {}
            for lipid_type in self._database_unique_resnames:
                occupancy_frames = []
                lipids_of_given_type = self._resids[self._resnames == lipid_type]
                for lipid_id in lipids_of_given_type:
                    values = self.contact_frames[res_id][lipid_id]
                    occupancy_frames.extend([fr for fr in values if fr not in occupancy_frames])
                occupancies_dict[res_id][lipid_type] = len(occupancy_frames) * 100 / self._nframes
        return occupancies_dict


class FittingFunctionMeta(type):
    """Metaclass for fitting functions."""

    def __init__(cls, name, bases, dct):
        if not hasattr(cls, "registry"):
            cls.registry = {}
        else:
            cls.registry[cls.name] = cls
        super().__init__(name, bases, dct)


class FittingFunction(metaclass=FittingFunctionMeta):
    """Base class for fitting functions."""

    name = None
    p0 = [1, 1, 1, 1]
    maxfev = 1000000

    def compute(self, x, *params):
        raise NotImplementedError("Subclasses must implement this method")

    def get_koff(self, popt):
        raise NotImplementedError("Subclasses must implement this method")

    def fit(self, x_data, y_data, **kwargs):
        if "p0" not in kwargs:
            kwargs["p0"] = self.p0
        if "maxfev" not in kwargs:
            kwargs["maxfev"] = self.maxfev
        popt, _ = curve_fit(self.compute, x_data, y_data, **kwargs)
        return popt
