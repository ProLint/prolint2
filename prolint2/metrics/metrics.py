r""":mod:`prolint2.metrics.metrics`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from typing import Callable, Iterable
import numpy as np

from prolint2.metrics.base import BaseMetric, Metric
from prolint2.metrics.registries import MetricRegistry
from prolint2.metrics.formatters import (
    DefaultOutputFormat,
    SingleOutputFormat,
    ProLintDashboardOutputFormat,
    CustomOutputFormat,
)


class UserDefinedMetric(BaseMetric):
    """
    A user-defined metric class.

    Parameters
    ----------
    custom_function : Callable[[Iterable[int]], float]
        A custom function to compute the metric.

    Attributes
    ----------
    name : str
        The name of the metric (set to 'custom').

    Methods
    -------
    compute_metric(contact_array)
        Compute the custom metric for a given contact array.

    Examples
    --------
    >>> custom_metric = UserDefinedMetric(my_custom_function)
    >>> result = custom_metric.compute_metric(contact_array)
    """

    name: str = "custom"

    def __init__(self, custom_function: Callable[[Iterable[int]], float]):
        super().__init__()
        self.custom_function = custom_function

    def compute_metric(self, contact_array):
        """
        Compute the custom metric for a given contact array.

        Parameters
        ----------
        contact_array : Iterable[int]
            An iterable representing contact data.

        Returns
        -------
        float
            The computed custom metric.

        Examples
        --------
        >>> result = custom_metric.compute_metric(contact_array)
        """
        return self.custom_function(contact_array)


class MeanMetric(BaseMetric):
    """
    A metric to compute the mean of contact values.

    Attributes
    ----------
    name : str
        The name of the metric (set to 'mean').

    Methods
    -------
    compute_metric(contact_array)
        Compute the mean of contact values for a given contact array.

    Examples
    --------
    >>> mean_metric = MeanMetric()
    >>> result = mean_metric.compute_metric(contact_array)
    """

    name: str = "mean"

    def compute_metric(self, contact_array):
        """
        Compute the mean of contact values for a given contact array.

        Parameters
        ----------
        contact_array : Iterable[int]
            An iterable representing contact data.

        Returns
        -------
        float
            The computed mean value.

        Examples
        --------
        >>> result = mean_metric.compute_metric(contact_array)
        """
        return np.mean(contact_array)


class SumMetric(BaseMetric):
    """
    A metric to compute the sum of contact values.

    Attributes
    ----------
    name : str
        The name of the metric (set to 'sum').

    Methods
    -------
    compute_metric(contact_array)
        Compute the sum of contact values for a given contact array.

    Examples
    --------
    >>> sum_metric = SumMetric()
    >>> result = sum_metric.compute_metric(contact_array)
    """

    name: str = "sum"

    def compute_metric(self, contact_array):
        """
        Compute the sum of contact values for a given contact array.

        Parameters
        ----------
        contact_array : Iterable[int]
            An iterable representing contact data.

        Returns
        -------
        int
            The computed sum of contact values.

        Examples
        --------
        >>> result = sum_metric.compute_metric(contact_array)
        """
        return np.sum(contact_array)


class MaxMetric(BaseMetric):
    """
    A metric to compute the maximum contact value.

    Attributes
    ----------
    name : str
        The name of the metric (set to 'max').

    Methods
    -------
    compute_metric(contact_array)
        Compute the maximum contact value for a given contact array.

    Examples
    --------
    >>> max_metric = MaxMetric()
    >>> result = max_metric.compute_metric(contact_array)
    """

    name: str = "max"

    def compute_metric(self, contact_array):
        """
        Compute the maximum contact value for a given contact array.

        Parameters
        ----------
        contact_array : Iterable[int]
            An iterable representing contact data.

        Returns
        -------
        int
            The maximum contact value.

        Examples
        --------
        >>> result = max_metric.compute_metric(contact_array)
        """
        return np.max(contact_array)


def create_metric(
    contacts,
    metrics,
    output_format=None,
    custom_function: Callable = None,
    metric_registry: MetricRegistry = None,
    lipid_type=None,
    **kwargs,
):
    """
    Create and configure a metric computation task.

    This function creates a metric computation task based on user-defined parameters. It returns a Metric object that can be used to compute metrics on contact data.

    Parameters
    ----------
    contacts : ContactData
        The contact data to analyze.
    metrics : list
        A list of metrics to compute, including 'max', 'sum', 'mean', or 'custom'.
    output_format : str, optional
        The output format for the computed metrics. Supported values include 'default', 'custom', 'single', or 'dashboard'.
    custom_function : Callable, optional
        A custom function to use when the 'custom' metric is selected.
    metric_registry : MetricRegistry, optional
        A registry containing metric classes.
    lipid_type : str, optional
        The type of lipid to consider in the metric computation.
    **kwargs
        Additional keyword arguments for output format classes.

    Returns
    -------
    Metric
        A Metric object configured with the specified parameters.

    Raises
    ------
    ValueError
        If an invalid output format is specified or if 'single' format is used with more than one metric.

    Examples
    --------
    >>> metrics = ['max', 'mean']
    >>> metric = create_metric(contacts, metrics, output_format='default')
    >>> metric.run()
    >>> results = metric.compute('max')
    >>> results = metric.compute('sum', lipid_type='DOPC')
    >>> metric.compute('median')  # Raises ValueError. Use `apply_function` instead.
    """
    if metric_registry is None:
        raise ValueError("A MetricRegistry instance must be provided.")

    output_format_classes = {
        "default": DefaultOutputFormat,
        "custom": CustomOutputFormat,
        "single": SingleOutputFormat,
        "dashboard": ProLintDashboardOutputFormat,
    }

    if len(metrics) != 1 and output_format == "single":
        raise ValueError(
            "The 'single' output format can only be used with a single metric."
        )

    if len(metrics) == 1 and output_format is None or output_format == "single":
        output_format_class = SingleOutputFormat()
    else:
        if output_format is None:
            output_format = "default"
        if output_format not in output_format_classes:
            raise ValueError(
                f"Invalid output format '{output_format}'. Supported output formats are {list(output_format_classes.keys())}"
            )

        output_format_class = output_format_classes[output_format](**kwargs)

    metric_objects = []
    for metric in metrics:
        metric_class = metric_registry.get_metric(metric)

        if metric == "custom":
            if custom_function is not None:
                metric_objects.append(metric_class(custom_function))
            else:
                raise ValueError(
                    "A custom function must be provided when using the 'custom' metric."
                )
        else:
            metric_objects.append(metric_class())

    return Metric(contacts, metric_objects, output_format_class, lipid_type)
