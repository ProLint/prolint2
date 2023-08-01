from typing import Callable, Iterable
import numpy as np

from prolint2.metrics.base import BaseMetric, Metric
from prolint2.metrics.registries import MetricRegistry
from prolint2.metrics.formatters import DefaultOutputFormat, SingleOutputFormat, ProLintDashboardOutputFormat, CustomOutputFormat

    
class UserDefinedMetric(BaseMetric):
    name: str = 'custom'
    def __init__(self, custom_function: Callable[[Iterable[int]], float]):
        super().__init__()
        self.custom_function = custom_function

    def compute_metric(self, contact_array):
        return self.custom_function(contact_array)

class MeanMetric(BaseMetric):
    name: str = 'mean'
    def compute_metric(self, contact_array):
        return np.mean(contact_array)

class SumMetric(BaseMetric):
    name: str = 'sum'
    def compute_metric(self, contact_array):
        return np.sum(contact_array)

class MaxMetric(BaseMetric):
    name: str = 'max'
    def compute_metric(self, contact_array):
        return np.max(contact_array)

def create_metric(contacts, metrics, output_format=None, custom_function: Callable = None, metric_registry: MetricRegistry=None, lipid_type=None, **kwargs):
    if metric_registry is None:
        raise ValueError("A MetricRegistry instance must be provided.")

    output_format_classes = {
        'default': DefaultOutputFormat,
        'custom': CustomOutputFormat,
        'single': SingleOutputFormat,
        'dashboard': ProLintDashboardOutputFormat,
    }

    if len(metrics) != 1 and output_format == 'single':
        raise ValueError("The 'single' output format can only be used with a single metric.")
    
    if len(metrics) == 1 and output_format is None or output_format == 'single':
        output_format_class = SingleOutputFormat()
    else:
        if output_format is None:
            output_format = 'default'
        if output_format not in output_format_classes:
            raise ValueError(f"Invalid output format '{output_format}'. Supported output formats are {list(output_format_classes.keys())}")

        output_format_class = output_format_classes[output_format](**kwargs)

    metric_objects = []
    for metric in metrics:
        metric_class = metric_registry.get_metric(metric)

        if metric == 'custom':
            if custom_function is not None:
                metric_objects.append(metric_class(custom_function))
            else:
                raise ValueError("A custom function must be provided when using the 'custom' metric.")
        else:
            metric_objects.append(metric_class())

    return Metric(contacts, metric_objects, output_format_class, lipid_type)
