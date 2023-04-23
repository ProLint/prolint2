from abc import ABC, abstractmethod
from prolint2.utils.formatters import OutputFormat, DefaultOutputFormat

class BaseMetric(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def compute_metric(self, contact_array):
        pass

class Metric(ABC):
    def __init__(self, contacts, metrics, output_format: OutputFormat = DefaultOutputFormat(), lipid_type=None, clear=True):
        self.contact_input = contacts
        if not isinstance(metrics, list):
            metrics = [metrics]
        self.metrics = metrics
        if clear:
            output_format.clear()
        self.output_format = output_format
        self.lipid_type = lipid_type

    def compute(self, dt=1, totaltime=1):
        multiplier = dt / totaltime
        for residue_id, lipid_dict in self.contact_input.items():
            for lipid_id, lipid_contacts in lipid_dict.items():
                if self.lipid_type is not None and self.lipid_type != lipid_id:
                    continue
                contact_array = list(lipid_contacts.values())

                if contact_array:
                    for metric in self.metrics:
                        value = metric.compute_metric(contact_array) * multiplier
                        self.output_format.store_result(residue_id, lipid_id, metric.__class__.__name__, value)
                else:
                    for metric in self.metrics:
                        self.output_format.store_result(residue_id, lipid_id, metric.__class__.__name__, 0)

        return self.output_format.get_result()
