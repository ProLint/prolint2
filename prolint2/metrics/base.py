from abc import ABC, abstractmethod
from scipy.optimize import curve_fit
from prolint2.metrics.formatters import OutputFormat, DefaultOutputFormat

from typing import Type
MetricRegistry = Type["registries.MetricRegistry"]

class BaseMetric(ABC):
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
    def __init__(self, contacts, metrics, output_format: OutputFormat = DefaultOutputFormat(), lipid_type=None, clear=True):
        self.contact_input = dict(sorted(contacts.contacts.items()))

        if not isinstance(metrics, list):
            metrics = [metrics]
        self.metrics = metrics
        if clear:
            output_format.clear()
        self.output_format = output_format
        self.lipid_type = lipid_type

    def compute(self, dt=1, totaltime=1):
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
                        self.output_format.store_result(residue_id, lipid_name, metric.__class__.__name__, value)
                else:
                    for metric in self.metrics:
                        self.output_format.store_result(residue_id, lipid_name, metric.__class__.__name__, 0)

        return self.output_format.get_result()

class FittingFunctionMeta(type):
    def __init__(cls, name, bases, dct):
        if not hasattr(cls, 'registry'):
            cls.registry = {}
        else:
            cls.registry[cls.name] = cls
        super().__init__(name, bases, dct)

class FittingFunction(metaclass=FittingFunctionMeta):
    name = None
    p0 = [1, 1, 1, 1]
    maxfev = 1000000

    def compute(self, x, *params):
        raise NotImplementedError("Subclasses must implement this method")

    def get_koff(self, popt):
        raise NotImplementedError("Subclasses must implement this method")

    def fit(self, x_data, y_data, **kwargs):
        if 'p0' not in kwargs:
            kwargs['p0'] = self.p0
        if 'maxfev' not in kwargs:
            kwargs['maxfev'] = self.maxfev
        popt, _ = curve_fit(self.compute, x_data, y_data, **kwargs)
        return popt
