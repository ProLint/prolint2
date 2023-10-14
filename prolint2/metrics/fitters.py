import numpy as np
from scipy.optimize import curve_fit
from prolint2.metrics.base import FittingFunction


class BiExpoFittingFunction(FittingFunction):
    name = "bi_expo"
    p0 = [1, 1.0, 1.0, 1.0]
    maxfev = 1000000

    def compute(self, x, k1, k2, A, B):
        exp1 = np.exp(np.clip(-k1 * x, None, 700))
        exp2 = np.exp(np.clip(-k2 * x, None, 700))
        return A * exp1 + B * exp2

    def get_koff(self, popt):
        ks = [abs(k) for k in popt[:2]]
        return np.min(ks)


class MonoExpoFittingFunction(FittingFunction):
    name = "mono_expo"
    p0 = [1, 1]
    maxfev = 1000000

    def compute(self, x, k, A):
        exp = np.exp(np.clip(-k * x, None, 700))
        return A * exp

    def get_koff(self, popt):
        return abs(popt[0])


class PolynomialFittingFunction(FittingFunction):
    name = "poly"
    p0 = [1, 1, 1, 1]
    maxfev = 1000000

    def __init__(self, degree=None):
        self.degree = degree if degree is not None else 1

    def compute(self, x, *params):
        return np.polyval(params, x)

    def get_koff(self, popt):
        # Calculate koff based on the polynomial coefficients
        raise NotImplementedError(
            "koff calculation for polynomial fitting function is not implemented yet"
        )

    def fit(self, x_data, y_data, **kwargs):
        degree = kwargs.pop("degree", None)
        if degree is not None:
            self.degree = degree
        if "p0" not in kwargs and self.degree is not None:
            kwargs["p0"] = [1] * (self.degree + 1)
        kwargs.pop("degree", None)
        popt, _ = curve_fit(self.compute, x_data, y_data, **kwargs)
        return popt


class FittingFunctionFactory:
    @staticmethod
    def get_fitting_function(name):
        try:
            return FittingFunction.registry[name]()
        except KeyError:
            raise ValueError(f"Invalid fitting function name: {name}")
