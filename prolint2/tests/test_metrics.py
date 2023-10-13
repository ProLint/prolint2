import numpy as np
import pytest
from prolint2 import Universe
from prolint2.metrics.aprox_contacts import AproxContacts
from prolint2.metrics.exact_contacts import ExactContacts
from prolint2.sampledata import GIRKDataSample

GIRK = GIRKDataSample()

from prolint2.metrics.base import BaseMetric, Metric
from prolint2.metrics.registries import MetricRegistry
from prolint2.metrics.formatters import (
    DefaultOutputFormat,
    SingleOutputFormat,
)
from prolint2.metrics.metrics import (
    UserDefinedMetric,
    create_metric,
)

from prolint2.metrics.base import FittingFunction
from prolint2.metrics.fitters import (
    FittingFunctionFactory,
    BiExpoFittingFunction,
    MonoExpoFittingFunction,
    PolynomialFittingFunction,
)
from prolint2.metrics.restime import SurvivalFunction, KoffCalculator


@pytest.fixture(scope="module")
def contact_store():
    # Create a fixture for the BaseContactStore class
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts = u.compute_contacts(cutoff=7)
    return AproxContacts(u, contacts.contact_frames)


@pytest.fixture(scope="module")
def contact_store_exact():
    # Create a fixture for the BaseContactStore class
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts = u.compute_contacts(cutoff=7)
    return ExactContacts(u, contacts.contact_frames)


class TestAproxContacts:
    def test_run(self, contact_store):
        # Test the run() method

        # Define test data
        lipid_resnames = ["POPS"]

        # Call the run() method
        contact_store.run(lipid_resnames)

        # Perform assertions
        assert contact_store._contacts  # Check that contacts have been stored

    def test_pooled_results(self, contact_store):
        # Test the pooled_results() method

        # Define test data
        target_lipid_name = "POPS"

        # Call the pooled_results() method
        results = contact_store.pooled_results(target_lipid_name)

        # Perform assertions
        assert results  # Check that results are not empty
        assert isinstance(results, dict)  # Check that results are of type dict

    def test_compute_metric(self, contact_store):
        # Test the compute_metric() method

        # Define test data
        metric = "mean"
        target_lipid_name = "POPS"

        # Call the compute_metric() method
        results = contact_store.compute_metric(metric, target_lipid_name)

        # Perform assertions
        assert results  # Check that results are not empty
        assert isinstance(results, dict)  # Check that results are of type dict

    def test_apply_function(self, contact_store):
        # Test the apply_function() method

        # Define test data
        target_lipid_name = "POPS"

        # Define a sample function to apply
        def sample_function(array):
            return np.mean(array)

        # Call the apply_function() method
        results = contact_store.apply_function(sample_function, target_lipid_name)

        # Perform assertions
        assert results  # Check that results are not empty
        assert isinstance(results, dict)  # Check that results are of type dict


class TestExactContacts:
    def test_run(self, contact_store_exact):
        # Test the run() method

        # Define test data
        lipid_resnames = ["POPS"]

        # Call the run() method
        contact_store_exact.run(lipid_resnames)

        # Perform assertions
        assert contact_store_exact._contacts

    def test_pooled_results(self, contact_store_exact):
        # Test the pooled_results() method

        # Define test data
        target_lipid_name = "POPS"

        # Call the pooled_results() method
        results = contact_store_exact.pooled_results(target_lipid_name)

        # Perform assertions
        assert results

    def test_compute_metric(self, contact_store_exact):
        # Test the compute_metric() method

        # Define test data
        metric = "mean"
        target_lipid_name = "POPS"

        # Call the compute_metric() method
        results = contact_store_exact.compute_metric(metric, target_lipid_name)

        # Perform assertions
        assert results

    def test_apply_function(self, contact_store_exact):
        # Test the apply_function() method

        # Define test data
        target_lipid_name = "POPS"

        # Define a sample function to apply
        def sample_function(array):
            return np.mean(array)

        # Call the apply_function() method
        results = contact_store_exact.apply_function(sample_function, target_lipid_name)

        # Perform assertions
        assert results

    def test_compute_lipid_durations(self, contact_store_exact):
        # Test the compute_lipid_durations() method

        # Define test data
        lipid_resname = "CHOL"
        contact_frame = contact_store_exact.contact_frames[
            41
        ]  # frame with cholesterol contacts

        # Call the compute_lipid_durations() method
        durations = contact_store_exact.compute_lipid_durations(
            contact_frame, lipid_resname
        )

        # Perform assertions
        assert durations


@pytest.fixture
def metric_registry():
    return MetricRegistry()  # You can customize this fixture based on your needs


@pytest.fixture
def contact_array():
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts_out = u.compute_contacts(cutoff=7)
    return contacts_out  # Example contact array for testing


@pytest.fixture
def custom_function():
    def custom_metric_function(contact_array):
        # Replace this with your custom metric calculation logic
        return sum(contact_array) / len(contact_array)

    return custom_metric_function


def test_create_metric_with_single_output_format(metric_registry, contact_array):
    metrics = ["mean"]
    output_format = "single"
    metric = create_metric(
        contact_array, metrics, output_format, metric_registry=metric_registry
    )
    assert isinstance(metric, Metric)
    assert isinstance(metric.output_format, SingleOutputFormat)


def test_create_metric_with_custom_metric(
    metric_registry, contact_array, custom_function
):
    metrics = ["custom"]
    output_format = "default"
    metric = create_metric(
        contact_array,
        metrics,
        output_format,
        custom_function,
        metric_registry=metric_registry,
    )
    assert isinstance(metric, Metric)
    assert isinstance(metric.output_format, DefaultOutputFormat)
    assert isinstance(metric.metrics[0], UserDefinedMetric)


def test_create_metric_without_custom_function_and_custom_metric(
    metric_registry, contact_array
):
    metrics = ["custom"]
    output_format = "default"
    with pytest.raises(ValueError):
        create_metric(
            contact_array, metrics, output_format, metric_registry=metric_registry
        )


def test_create_metric_with_invalid_output_format(metric_registry, contact_array):
    metrics = ["mean", "max"]
    output_format = "invalid_format"
    with pytest.raises(ValueError):
        create_metric(
            contact_array, metrics, output_format, metric_registry=metric_registry
        )


def test_create_metric_with_multiple_metrics(metric_registry, contact_array):
    metrics = ["mean", "sum", "max"]
    output_format = "default"
    metric = create_metric(
        contact_array, metrics, output_format, metric_registry=metric_registry
    )
    assert isinstance(metric, Metric)
    assert isinstance(metric.output_format, DefaultOutputFormat)
    assert len(metric.metrics) == 3
    assert all(isinstance(metric_obj, BaseMetric) for metric_obj in metric.metrics)


class TestFittingFunctions:
    def test_bi_expo_compute(self):
        fitting_function = BiExpoFittingFunction()
        x = np.array([1, 2, 3])
        k1, k2, A, B = 1, 2, 3, 4
        expected_result = A * np.exp(np.clip(-k1 * x, None, 700)) + B * np.exp(
            np.clip(-k2 * x, None, 700)
        )
        assert np.allclose(fitting_function.compute(x, k1, k2, A, B), expected_result)

    def test_bi_expo_get_koff(self):
        fitting_function = BiExpoFittingFunction()
        popt = np.array([1, 2, 3, 4])
        expected_result = np.min(np.abs(popt[:2]))
        assert fitting_function.get_koff(popt) == expected_result

    def test_mono_expo_compute(self):
        fitting_function = MonoExpoFittingFunction()
        x = np.array([1, 2, 3])
        k, A = 1, 2
        expected_result = A * np.exp(np.clip(-k * x, None, 700))
        assert np.allclose(fitting_function.compute(x, k, A), expected_result)

    def test_mono_expo_get_koff(self):
        fitting_function = MonoExpoFittingFunction()
        popt = np.array([1, 2])
        expected_result = np.abs(popt[0])
        assert fitting_function.get_koff(popt) == expected_result

    def test_poly_compute(self):
        fitting_function = PolynomialFittingFunction()
        x = np.array([1, 2, 3])
        params = [1, 2, 3]
        expected_result = np.polyval(params, x)
        assert np.allclose(fitting_function.compute(x, *params), expected_result)

    def test_poly_get_koff(self):
        fitting_function = PolynomialFittingFunction()
        popt = np.array([1, 2, 3])
        with pytest.raises(NotImplementedError):
            fitting_function.get_koff(popt)

    def test_fitting_function_factory_get_fitting_function_valid_name(self):
        fitting_function = FittingFunctionFactory.get_fitting_function("bi_expo")
        assert isinstance(fitting_function, BiExpoFittingFunction)

    def test_fitting_function_factory_get_fitting_function_invalid_name(self):
        with pytest.raises(ValueError):
            FittingFunctionFactory.get_fitting_function("invalid_name")


class TestSurvivalFunction:
    @pytest.fixture
    def survival_function(self):
        durations = [1, 2, 3, 4, 5]
        t_total = 10
        delta_t_list = [0, 1, 2, 3, 4, 5]
        return SurvivalFunction(durations, t_total, delta_t_list)

    def test_calc_survival_value(self, survival_function):
        delta_t = 2
        expected_result = 0.5
        assert survival_function._calc_survival_value(delta_t) == pytest.approx(
            expected_result, 0.1
        )


class TestKoffCalculator:
    @pytest.fixture
    def koff_calculator(self):
        durations = [0.39, 0.79, 0.39, 0.39]
        t_total = 4.8
        timestep = 0.4
        fitting_func_name = "bi_expo"
        return KoffCalculator(
            durations, t_total, timestep, fitting_func_name=fitting_func_name
        )

    def test_is_empty_or_zeros(self, koff_calculator):
        assert not koff_calculator._is_empty_or_zeros([1, 2, 3])
        assert koff_calculator._is_empty_or_zeros([])
        assert koff_calculator._is_empty_or_zeros(np.array([0, 0, 0]))

    def test_calculate_koff(self, koff_calculator):
        expected_res_time = 0.384892
        expected_koff = 2.59812
        koff_calculator.calculate_koff()
        assert koff_calculator.res_time == pytest.approx(expected_res_time, 0.1)
        assert koff_calculator.koff == pytest.approx(expected_koff, 0.1)

    def test_invalid_fitting_func_name(self):
        with pytest.raises(ValueError):
            durations = [1, 2, 3, 4, 5]
            t_total = 10
            timestep = 1
            fitting_func_name = "invalid_func"
            KoffCalculator(
                durations, t_total, timestep, fitting_func_name=fitting_func_name
            )

    def test_calculate_koff_empty_durations(self):
        durations = []
        t_total = 10
        timestep = 1
        fitting_func_name = "bi_expo"
        koff_calculator = KoffCalculator(
            durations, t_total, timestep, fitting_func_name=fitting_func_name
        )
        assert koff_calculator.res_time == 0
        assert koff_calculator.koff == 0


class TestFittingFunctionFactory:
    def test_get_fitting_function(self):
        fitting_func_name = "bi_expo"
        fitting_func = FittingFunctionFactory.get_fitting_function(fitting_func_name)
        assert isinstance(fitting_func, FittingFunction)
        assert fitting_func_name in FittingFunction.registry
