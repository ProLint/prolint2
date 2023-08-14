import numpy as np
from collections import defaultdict


def use_1d_script_template(code_body, tail):
    template = """
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from prolint2 import Universe
from prolint2.plot.plot_1d import Plotter1D
from prolint2.plot.utils import *

## seaborn config for paper quality plots
import seaborn as sns
sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")

{} 

if __name__ == "__main__":
    # Define your data structure here: including the Prolint2 universe, the contacts and the metrics as needed
    # Example:
    from prolint2.sampledata import GIRKDataSample
    from prolint2.metrics.metrics import Metric, MeanMetric, SumMetric, MaxMetric
    GIRK = GIRKDataSample()
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts = u.compute_contacts(cutoff=7)

{}
    """.format(
        code_body, tail
    )
    return template


def get_metric_list_by_residues(universe, metric, lipid=None, metric_name=None):
    """Get a list of metric values for a given lipid and metric name"""
    resids = universe.query.residues.resids
    metric_list = []
    for resi in resids:
        if resi in metric.keys():
            if metric_name:
                if lipid:
                    try:
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        metric_list.append(0)
                else:
                    # get the lipid names available
                    lipids_dict = {
                        i: lipid
                        for i, lipid in enumerate(universe.database.unique_resnames)
                    }
                    lipid = int(
                        input(
                            "Please enter the lipid name that you want to analyze: \n {}.".format(
                                lipids_dict
                            )
                        )
                    )
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available".format(lipid))
                    if metric_name not in metric[resi][lipid].keys():
                        # break the loops and print error message
                        print(
                            "The metric {} is not available for the lipid {}.".format(
                                metric_name, lipid
                            )
                        )
                        break
                    try:
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        print(
                            "The metric {} is not available for the lipid {}.".format(
                                metric_name, lipid
                            )
                        )
            else:
                if lipid:
                    if lipid not in universe.database.unique_resnames:
                        # break the loops and print error message
                        print("The lipid {} is not available.".format(lipid))
                        break
                else:
                    # get the lipid names available
                    lipids_dict = {
                        i: lipid
                        for i, lipid in enumerate(universe.database.unique_resnames)
                    }
                    lipid = int(
                        input(
                            "Please enter the lipid name that you want to analyze: \n {}.".format(
                                lipids_dict
                            )
                        )
                    )
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available.".format(lipid))
                # get the metric names available
                metric_names = list(metric[resi][lipid].keys())
                inc = 1
                while len(metric_names) == 0:
                    metric_names = list(metric[resi + inc][lipid].keys())
                    inc += 1
                metric_names_dict = {
                    i: metric_name for i, metric_name in enumerate(metric_names)
                }
                metric_name = int(
                    input(
                        "Please enter the metric name that you want to analyze: \n {}.".format(
                            metric_names_dict
                        )
                    )
                )
                try:
                    metric_name = metric_names_dict[metric_name]
                except KeyError:
                    print("The metric {} is not available.".format(metric_name))
                try:
                    metric_list.append(metric[resi][lipid][metric_name])
                except KeyError:
                    metric_list.append(0)
        else:
            metric_list.append(0)
    return np.array(metric_list)


class AxisIndex:
    """Build axes for logo figure."""

    def __init__(self, residue_index, logos, interactions, length, gap):
        self.page_idx = 0
        self.length = length
        self.gap = gap
        self.residue_index = residue_index
        self.logos = logos
        self.interactions = interactions
        self.axis_start = (residue_index[0] // length) * length
        self.breaks = defaultdict(list)
        self.breaks[self.page_idx].append([])
        self.gray_areas = defaultdict(list)

    def fill_missing(self, start_value, end_value):
        for xloci in np.arange(start_value, end_value + 1):
            self.breaks[self.page_idx][-1].append((xloci, "A", 0))
        self.gray_areas[self.page_idx].append(
            (len(self.breaks[self.page_idx]) - 1, start_value, end_value)
        )

    def new_axis(self, pointer):
        self.breaks[self.page_idx].append([])
        self.axis_start = self.residue_index[pointer]
        self.breaks[self.page_idx][-1].append(
            (
                self.residue_index[pointer],
                self.logos[pointer],
                self.interactions[pointer],
            )
        )

    def new_page(self, pointer):
        if len(self.breaks[self.page_idx][-1]) < self.length:
            self.fill_missing(
                self.axis_start + len(self.breaks[self.page_idx][-1]),
                self.axis_start + self.length - 1,
            )
        self.page_idx += 1
        self.breaks[self.page_idx].append([])
        self.axis_start = (self.residue_index[pointer] // self.length) * self.length
        if self.axis_start != self.residue_index[pointer]:
            self.fill_missing(self.axis_start, self.residue_index[pointer] - 1)
        self.breaks[self.page_idx][-1].append(
            (
                self.residue_index[pointer],
                self.logos[pointer],
                self.interactions[pointer],
            )
        )

    def new_gap(self, pointer):
        gray_start = self.residue_index[pointer - 1] + 1
        for xloci in np.arange(
            self.residue_index[pointer - 1] + 1, self.residue_index[pointer]
        ):
            if xloci - self.axis_start < self.length:
                self.breaks[self.page_idx][-1].append((xloci, "A", 0))
            else:
                self.gray_areas[self.page_idx].append(
                    (len(self.breaks[self.page_idx]) - 1, gray_start, xloci - 1)
                )
                self.breaks[self.page_idx].append([])
                self.breaks[self.page_idx][-1].append((xloci, "A", 0))
                self.axis_start = xloci
                gray_start = xloci
        self.gray_areas[self.page_idx].append(
            (
                len(self.breaks[self.page_idx]) - 1,
                gray_start,
                self.residue_index[pointer] - 1,
            )
        )
        self.breaks[self.page_idx][-1].append(
            (
                self.residue_index[pointer],
                self.logos[pointer],
                self.interactions[pointer],
            )
        )

    def sort(self):
        end = False
        if self.axis_start != self.residue_index[0]:
            self.fill_missing(self.axis_start, self.residue_index[0] - 1)
        self.breaks[self.page_idx][-1].append(
            (self.residue_index[0], self.logos[0], self.interactions[0])
        )
        pointer = 1
        while not end:
            if (
                self.residue_index[pointer] - self.residue_index[pointer - 1] == 1
                and self.residue_index[pointer] - self.axis_start < self.length
            ):
                self.breaks[self.page_idx][-1].append(
                    (
                        self.residue_index[pointer],
                        self.logos[pointer],
                        self.interactions[pointer],
                    )
                )
                pointer += 1
            elif (
                self.residue_index[pointer] - self.residue_index[pointer - 1] == 1
                and self.residue_index[pointer] - self.axis_start >= self.length
            ):
                self.new_axis(pointer)
                pointer += 1
            elif self.residue_index[pointer] - self.residue_index[pointer - 1] < 0:
                self.new_page(pointer)
                pointer += 1
            elif (
                1
                < self.residue_index[pointer] - self.residue_index[pointer - 1]
                <= self.gap
            ):
                self.new_gap(pointer)
                pointer += 1
            elif (
                self.residue_index[pointer] - self.residue_index[pointer - 1] > self.gap
            ):
                self.new_page(pointer)
                pointer += 1
            if pointer == len(self.residue_index):
                end = True
        if len(self.breaks[self.page_idx][-1]) < self.length:
            self.fill_missing(
                self.axis_start + len(self.breaks[self.page_idx][-1]),
                self.axis_start + self.length - 1,
            )

def get_metrics_for_radar(metrics, metric_names, resIDs=[], lipid=None):
    metrics_radar = {}
    for resi in resIDs:
        metrics_radar[resi] = []
        if resi in metrics.keys():
            for metric in metric_names:
                if metric in metrics[resi][lipid].keys():
                    metrics_radar[resi].append(metrics[resi][lipid][metric])
                else:
                    metrics_radar[resi].append(0)
        else:
            for metric in metric_names:
                metrics_radar[resi].append(0)

    return metrics_radar, metric_names