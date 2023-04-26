from typing import List, Type
import importlib
import inspect

from prolint2.metrics.base import BaseMetric

class MetricRegistry:
    def __init__(self):
        self._metrics = {}

    def register(self, name: str, metric_class: Type[BaseMetric]):
        if name in self._metrics:
            raise ValueError(f"A metric with the name '{name}' is already registered.")
        self._metrics[name] = metric_class

    def get_metric(self, name: str) -> Type[BaseMetric]:
        if name not in self._metrics:
            raise ValueError(f"No metric found with the name '{name}'.")
        return self._metrics[name]

    def get_registered_names(self) -> List[str]:
        return list(self._metrics.keys())
    

def auto_register_metrics(metric_registry: MetricRegistry, module_name: str):
    module = importlib.import_module(module_name)
    for _, obj in inspect.getmembers(module):
        if inspect.isclass(obj) and issubclass(obj, BaseMetric) and obj != BaseMetric:
            print ('Registering metric: ', obj.name, ' from module: ', module_name, '...')
            metric_name = obj.name
            metric_registry.register(metric_name, obj)
