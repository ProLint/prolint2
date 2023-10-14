Tutorial #3: Introduction to the **Contacts** and **Metrics** objects
=====================================================================

This tutorial will walk you through the essential concepts and functionalities for computing various metrics related to lipid-protein interactions in your molecular dynamics simulations.

Prerequisites
-------------

Before we begin, make sure you have **prolint2** installed in your Python environment. 

.. code-block:: python

    from typing import List, Iterable
    import numpy as np
    from prolint2 import Universe
    from prolint2.sampledata import GIRKDataSample
    GIRK = GIRKDataSample()
    ts = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts = ts.compute_contacts(cutoff=7)

Contact Output
--------------

The **Contacts** object provides you with a all the information about lipid-protein interactions in your system. By using the `contacts` object, you can access non-formatted contact data, which includes a nested structure containing all contact information. The key attributes are `contacts.contact_frames` and `contacts.contacts`. These attributes contain detailed information about the contacts formed during your trajectory.

.. code-block:: python

    # Access non-formatted contact output
    contacts.contact_frames
    contacts.contacts

Computing Contact Metrics
-------------------------

**prolint2** offers various metrics to analyze and quantify lipid-protein interactions. You can easily compute these metrics using the provided classes. Here are a few examples:

* Computing Mean, Sum, and Max Metrics

.. code-block:: python

    from prolint2.metrics.metrics import MeanMetric, SumMetric, MaxMetric

    # Create instances of the metrics you want to compute
    mean_instance = MeanMetric()
    sum_instance = SumMetric()
    max_instance = MaxMetric()

    # Compute the metrics using the `Metric` class
    mean_contacts = Metric(contacts, mean_instance).compute()
    sum_contacts = Metric(contacts, sum_instance).compute()
    max_contacts = Metric(contacts, max_instance).compute()

    # Access the metrics for a specific residue (e.g., residue 14)
    mean_contacts[14]

* Defining Custom Metrics

You can also define your own custom metrics by creating a class that inherits from the `BaseMetric` class. Here are two examples:

.. code-block:: python

    from prolint2.metrics.base import BaseMetric

    class ScaleAndMeanMetric(BaseMetric):
        name: str = 'scale'

        def compute_metric(self, contact_array: Iterable) -> float:
            return np.mean(contact_array) * 2

    class RandomWeightedMeanMetric(BaseMetric):
        name: str = 'weighted_mean'

        def compute_metric(self, contact_array: Iterable) -> float:
            return np.average(contact_array, weights=np.random.rand(len(contact_array)))

    # Create instances of the custom metrics and compute them
    scale_and_mean_instance = Metric(contacts, ScaleAndMeanMetric())
    scale_and_mean_contacts = scale_and_mean_instance.compute()

    weighted_mean_instance = Metric(contacts, RandomWeightedMeanMetric())
    weighted_mean_contacts = weighted_mean_instance.compute()

* Using User-Defined Metrics

You can also define a custom metric using your own function and then compute it. ProLint provides a `UserDefinedMetric` class for this purpose.

.. code-block:: python

    from prolint2.metrics.metrics import UserDefinedMetric

    # Define a custom metric function
    def custom_user_function(contact_array: Iterable) -> float:
        return np.mean(contact_array) * 10

    # Create an instance of the UserDefinedMetric class with your custom function
    user_metric_instance = UserDefinedMetric(custom_user_function)

    # Compute the user-defined metric
    user_metric = Metric(contacts, user_metric_instance)
    user_metric_contacts = user_metric.compute()

* Appending Metrics and Computing Multiple Metrics

You can choose to append results to the metric output or compute multiple metrics at once.

.. code-block:: python

    # Append results to the metric output
    metric_instance = Metric(contacts, MeanMetric())  # By default, `clear` is True, clearing any existing metrics
    contacts_out = metric_instance.compute()

    metric_instance = Metric(contacts, SumMetric(), clear=False)  # Set `clear` to False to keep existing metrics
    contacts_out = metric_instance.compute()

    # Compute multiple metrics at once
    metric_instances_list = [MeanMetric(), SumMetric(), MaxMetric()]
    metric_instance = Metric(contacts, metric_instances_list)  # By default, `clear` is True, clearing any existing metrics
    contacts_out = metric_instance.compute()

Output Formats
--------------

ProLint allows you to choose from various output formats, depending on your needs. The default output format is `DefaultOutputFormat`, but you can use other formats like `SingleOutputFormat`, `CustomOutputFormat`, and `ProLintDashboardOutputFormat`.

.. code-block:: python

    from prolint2.metrics.formatters import DefaultOutputFormat, SingleOutputFormat, CustomOutputFormat, ProLintDashboardOutputFormat

    # Example of using a custom output format
    metric_instances_list = [MeanMetric(), SumMetric(), MaxMetric()]
    metric_instance = Metric(contacts, metric_instances_list, output_format=CustomOutputFormat())
    contacts_out = metric_instance.compute()

    # The `ProLintDashboardOutputFormat` is used by the ProLint Dashboard and requires residue names and IDs
    input_dict = {
        'residue_names': ts.query.residues.resnames,
        'residue_ids': ts.query.residues.resids
    }

    metric_instance = Metric(contacts, MeanMetric(), output_format=ProLintDashboardOutputFormat(**input_dict))
    contacts_out = metric_instance.compute()

    # For single metrics, you can use the `SingleOutputFormat`
    metric_instance = Metric(contacts, MeanMetric(), output_format=SingleOutputFormat())
    contacts_out = metric_instance.compute()

Creating and Adding Metrics to the Registry
-------------------------------------------

You can add your custom metrics to the ProLint metric registry for easy access. Here's how you can do it:

.. code-block:: python

    # Add your custom metric class to the registry
    registry.register('scaled_mean', ScaleAndMeanMetric)

    # You can now use your custom metric by referring to it by name
    metric_instance = create_metric(contacts, metrics=['scaled_mean', 'max', 'mean'], metric_registry=registry, output_format='default')
    contacts_out = metric_instance.compute()

Converting Output Formats
-------------------------

You can convert between different output formats using the provided converters, such as `DefaultToSingleConverter` and `CustomToSingleConverter`.

.. code-block:: python
    from prolint2.metrics.converters import DefaultToSingleConverter, CustomToSingleConverter

    # Convert from the default output format to the single output format
    metric_instance = create_metric(contacts, metrics=['scaled_mean', 'max', 'mean'], metric_registry=registry, output_format='default')
    contacts_out = metric_instance.compute()

    # Extract a single metric from the converted output
    extract_single_metric = DefaultToSingleConverter(contacts_out, 'scaled_mean', registry).convert().get_result()

    # Convert from the custom output format to the single output format
    metric_instance = create_metric(contacts, metrics=['scaled_mean', 'max', 'mean'], metric_registry=registry, output_format='custom')
    contacts_out = metric_instance.compute()

    # Extract a single metric from the converted output
    extract_single_metric = CustomToSingleConverter(contacts_out, 0, registry).convert().get_result()

Convenience Function: `create_metric`
-------------------------------------

The `create_metric` function simplifies the process of creating a metric instance and computing the metric in one step. 

.. code-block:: python

    from prolint2.metrics.metrics import create_metric

    # Creating a metric with predefined metrics and output format
    registry = ts.registry
    metric_instance = create_metric(
        contacts,
        metrics=['mean', 'sum', 'max'],
        metric_registry=registry,
        output_format='default'
    )

    contacts_out = metric_instance.compute()

    # Creating a metric with a custom function
    def custom_function(contact_array: Iterable) -> float:
        return np.mean(contact_array) * 10

    metric_instance = create_metric(
        contacts,
        metrics=['custom'],
        custom_function=custom_function,
        metric_registry=registry,
        output_format='default'
    )

    contacts_out = metric_instance.compute()