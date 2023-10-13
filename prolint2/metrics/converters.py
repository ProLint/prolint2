from abc import ABC, abstractmethod
from typing import Union

from prolint2.metrics.formatters import OutputFormat, SingleOutputFormat
from prolint2.metrics.registries import MetricRegistry


class BaseOutputFormatConverter(ABC):
    def __init__(
        self,
        output_format: OutputFormat,
        metric_type: Union[str, int],
        metric_registry: MetricRegistry,
    ):
        self.output_format = output_format
        self.metric_type = metric_type
        self.metric_registry = metric_registry

    @abstractmethod
    def convert(self) -> OutputFormat:
        pass


class DefaultToSingleConverter(BaseOutputFormatConverter):
    def convert(self) -> SingleOutputFormat:
        single_output_format = SingleOutputFormat()
        metric_class_name = self.metric_registry.get_metric(self.metric_type).__name__

        if not isinstance(self.metric_type, str):
            raise ValueError(
                "The metric type must be an string when using the DefaultToSingleConverter"
            )

        for residue_id, lipid_data in self.output_format.items():
            for lipid_id, metric_data in lipid_data.items():
                single_output_format.store_result(
                    residue_id,
                    lipid_id,
                    self.metric_type,
                    metric_data[metric_class_name],
                )

        return single_output_format


class CustomToSingleConverter(BaseOutputFormatConverter):
    def convert(self) -> SingleOutputFormat:
        single_output_format = SingleOutputFormat()

        if not isinstance(self.metric_type, int):
            raise ValueError(
                "The metric type must be an integer when using the CustomToSingleConverter"
            )

        for residue_id, lipid_data in self.output_format.items():
            for lipid_id, metric_values in lipid_data.items():
                single_output_format.store_result(
                    residue_id,
                    lipid_id,
                    self.metric_type,
                    metric_values[self.metric_type],
                )

        return single_output_format
