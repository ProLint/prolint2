from typing import List, Type
import importlib
import inspect

import logging

from prolint2.metrics.base import BaseMetric


class MetricRegistry:
    def __init__(self):
        self._metrics = {}
        self.module_name = "prolint2.metrics.metrics"

        module = importlib.import_module(self.module_name)
        for _, obj in inspect.getmembers(module):
            if (
                inspect.isclass(obj)
                and issubclass(obj, BaseMetric)
                and obj != BaseMetric
            ):
                metric_name = obj.name
                self.register(metric_name, obj)

        # for metric_class in BaseMetric.__subclasses__():
        #     metric_class._register(self)

    def register(self, name: str, metric_class: Type[BaseMetric]):
        if name in self._metrics:
            logging.warning(
                lambda: "Metric with name '%s' already exists in registry.", name
            )
        self._metrics[name] = metric_class

    def get_metric(self, name: str) -> Type[BaseMetric]:
        if name not in self._metrics:
            raise ValueError(f"No metric found with the name '{name}'.")
        return self._metrics[name]

    def get_registered_names(self) -> List[str]:
        return list(self._metrics.keys())
