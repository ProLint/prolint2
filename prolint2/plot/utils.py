import numpy as np
from collections import defaultdict


def use_1d_script_template(code_body, tail):
    template = """
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from prolint2 import Universe
from prolint2.plot.plp import Plotter
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


def compute_density(universe, lipids, size_in_mb):
    frames = universe.trajectory.n_frames
    selection_string = "resname {}".format(lipids[0])
    if len(lipids) > 1:
        for l in lipids[1:]:
            selection_string += " or resname {}".format(l)

    lipids_ag = universe.select_atoms(selection_string)
    for ts in universe.trajectory:
        if ts.frame == 0:
            lipids_xyz = lipids_ag.positions
        else:
            lipids_xyz = np.append(lipids_xyz, lipids_ag.positions)

    # Convert size from MB to Bytes
    size = size_in_mb * (1024**2)  # 1 MB = 1024 * 1024 bytes

    def slice_array(arr, slice_by):
        mask = np.ones_like(arr, dtype=bool)
        print(mask.shape)
        mask[:: int(slice_by)] = False
        print(mask[:: int(slice_by)].shape[0])
        xshape = int(arr.shape[0] - mask[:: int(slice_by)].shape[0])
        print(xshape)
        arr = arr[mask].reshape(xshape, arr.shape[1], 3)
        return arr

    lipids_xyz = lipids_xyz.reshape((frames, lipids_ag.atoms.n_atoms, 3))

    while lipids_xyz.nbytes > size:
        lipids_xyz = slice_array(lipids_xyz, 10)

    lipids_xyz = lipids_xyz.reshape(lipids_xyz.shape[0] * lipids_xyz.shape[1], 3)

    return lipids_xyz
