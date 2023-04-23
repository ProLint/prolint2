
from collections import defaultdict
from abc import ABC, abstractmethod

class OutputFormat(ABC):
    def __init__(self):
        self.class_name = self.__class__.__name__
        # pass

    @abstractmethod
    def store_result(self, residue_id, lipid_id, metric_name, value):
        pass

    @abstractmethod
    def get_result(self):
        pass

class DefaultOutputFormat(OutputFormat):
    def __init__(self):
        super().__init__()
        self.results = defaultdict(lambda: defaultdict(dict))

    def store_result(self, residue_id, lipid_id, metric_name, value):
        self.results[residue_id][lipid_id][metric_name] = value

    def get_result(self):
        return self.results

class CustomOutputFormat(OutputFormat):
    def __init__(self):
        super().__init__()
        self.results = defaultdict(lambda: defaultdict(list))

    def store_result(self, residue_id, lipid_id, metric_name, value):
        self.results[residue_id][lipid_id].append(value)

    def get_result(self):
        return self.results

class SingleOutputFormat(OutputFormat):
    def __init__(self):
        super().__init__()
        self.results = defaultdict(dict)

    def store_result(self, residue_id, lipid_id, metric_name, value):
        self.results[residue_id][lipid_id] = value

    def get_result(self):
        return self.results

class ProLintDashboardOutputFormat(OutputFormat):
    def __init__(self, residue_names=None, residue_ids=None, **kwargs):
        super().__init__()
        self.results = defaultdict(list)
        self.residue_names = residue_names
        self.residue_ids = residue_ids

    def store_result(self, residue_id, lipid_id, metric_name, value):
        if not value > 0:
            return
        
        self.results[lipid_id].append(
            {
                "residue": f"{self.residue_names[residue_id]} {self.residue_ids[residue_id]}",
                "value": float("{:.2f}".format(value)),
            }
        )

    def get_result(self):
        return self.results
