from abc import ABC, abstractmethod
from MDAnalysis.analysis.base import AnalysisBase


class ContactComputerBase(AnalysisBase, ABC):
    # @abstractmethod
    # def _compute_pairs(self):
    #     pass

    # @abstractmethod
    # def _compute(self):
    #     pass

    def __add__(self, other):
        pass

    def intersection(self, other):
        pass

    def union(self, other):
        pass
