import numpy as np

from prolint2.metrics.base import FittingFunction
from prolint2.metrics.fitters import FittingFunctionFactory


class SurvivalFunction:
    def __init__(self, durations, t_total, delta_t_list):
        self.durations = durations
        self.t_total = t_total
        self.delta_t_list = delta_t_list
        self.num_of_contacts = len(durations)
        self.survival_func = self.calculate()

    def _calc_survival_value(self, delta_t):
        filtered_durations = [
            res_time for res_time in self.durations if res_time >= delta_t
        ]
        sum_res_time = sum(filtered_durations) - delta_t * len(filtered_durations)
        denominator = (self.t_total - delta_t) * self.num_of_contacts

        if delta_t != 0:
            denominator *= self.survival_func0

        return float(sum_res_time) / denominator if denominator != 0 else 0

    def calculate(self):
        survival_func = {}
        for delta_t in self.delta_t_list:
            survival_value = self._calc_survival_value(delta_t)

            if delta_t == 0:
                survival_func[delta_t] = 1
                self.survival_func0 = survival_value
            else:
                survival_func[delta_t] = survival_value

        return survival_func


class KoffCalculator:
    def __init__(
        self, durations, t_total, timestep, fitting_func_name="bi_expo", **kwargs
    ):
        self.durations = durations
        self.t_total = t_total
        self.timestep = timestep
        self.delta_t_list = np.arange(0, t_total, timestep)
        self.kwargs = kwargs

        if self._is_empty_or_zeros(self.durations):
            self.res_time, self.koff = 0, 0
            return

        if fitting_func_name not in FittingFunction.registry:
            func_names = ", ".join(FittingFunction.registry.keys())
            raise ValueError(
                f"Invalid fitting_func_name: {fitting_func_name}. Valid names are: {func_names}"
            )

        self.fitting_func = FittingFunctionFactory.get_fitting_function(
            fitting_func_name
        )
        self.survival_func = SurvivalFunction(
            self.durations, np.max(self.t_total), self.delta_t_list
        ).survival_func
        self.res_time, self.koff = self.calculate_koff()

    def _is_empty_or_zeros(self, array):
        return len(array) == 0 or np.all(array == 0)

    def calculate_koff(self):
        survival_rates = np.nan_to_num(
            [self.survival_func[delta_t] for delta_t in self.delta_t_list]
        )  # TODO: check if nan_to_num is needed
        popt = self.fitting_func.fit(
            np.array(self.delta_t_list), np.array(survival_rates), **self.kwargs
        )
        koff = self.fitting_func.get_koff(popt)
        res_time = 1 / koff
        return res_time, koff
