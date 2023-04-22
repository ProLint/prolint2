from collections import defaultdict
from typing import Callable, Iterable
import numpy as np
from abc import ABC, abstractmethod

class Metric(ABC):
    def __init__(self, ts):
        self.ts = ts

    def compute(self, apply_multiplier=True):
        multiplier = self.ts.dt / self.ts.totaltime if apply_multiplier else 1
        results = defaultdict(dict)
        for residue_id, lipid_dict in self.ts.contacts.items():
            for lipid_id, lipid_contacts in lipid_dict.items():
                contact_array = list(lipid_contacts.values())

                if contact_array:
                    results[residue_id][lipid_id] = self.compute_metric(contact_array) * multiplier
                else:
                    results[residue_id][lipid_id] = 0

        return results

    @abstractmethod
    def compute_metric(self, contact_array):
        pass

class UserDefinedMetric(Metric):
    def __init__(self, ts, custom_function: Callable[[Iterable[int]], float]):
        super().__init__(ts)
        self.custom_function = custom_function

    def compute_metric(self, contact_array):
        return self.custom_function(contact_array)

class MeanMetric(Metric):
    def compute_metric(self, contact_array):
        return np.mean(contact_array)

class SumMetric(Metric):
    def compute_metric(self, contact_array):
        return np.sum(contact_array)

class MaxMetric(Metric):
    def compute_metric(self, contact_array):
        return np.max(contact_array)

def create_metric(ts, metric, custom_function=None):
    metric_classes = {
        'mean': MeanMetric,
        'sum': SumMetric,
        'max': MaxMetric,
        'custom': UserDefinedMetric,
    }

    if metric not in metric_classes:
        raise ValueError(f"Invalid metric '{metric}'. Supported metrics are {list(metric_classes.keys())}")

    metric_class = metric_classes[metric]

    if metric == 'custom' and custom_function is not None:
        return metric_class(ts, custom_function)
    else:
        return metric_class(ts)
