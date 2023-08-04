import numpy as np


def use_1d_script_template(code_body, tail):
    template = """
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
from prolint2 import Universe
from prolint2.plot.plot_1d import Plotter1D
from prolint2.plot.utils import *

{} 

if __name__ == "__main__":
    # Define your data structure here: including the Prolint2 universe, the contacts and the metrics as needed
    # Example:
    from prolint2.sampledata import GIRKDataSample
    from prolint2.metrics.metrics import Metric, MeanMetric
    GIRK = GIRKDataSample()
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    contacts = u.compute_contacts(cutoff=7)
    mean_instance = MeanMetric()
    metric_instance = Metric(contacts, mean_instance)
    mean_contacts = metric_instance.compute()
    mean_contacts
{}
    """.format(code_body, tail)
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
                    lipids_dict = {i: lipid for i, lipid in enumerate(universe.database.unique_resnames)}
                    lipid = int(input("Please enter the lipid name that you want to analyze: \n {}.".format(lipids_dict)))
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available".format(lipid))
                    if metric_name not in metric[resi][lipid].keys():
                        # break the loops and print error message
                        print("The metric {} is not available for the lipid {}.".format(metric_name, lipid))
                        break
                    try:
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        print("The metric {} is not available for the lipid {}.".format(metric_name, lipid))
            else:
                if lipid:
                    if lipid not in universe.database.unique_resnames:
                        # break the loops and print error message
                        print("The lipid {} is not available.".format(lipid))
                        break
                else:
                    # get the lipid names available
                    lipids_dict = {i: lipid for i, lipid in enumerate(universe.database.unique_resnames)}
                    lipid = int(input("Please enter the lipid name that you want to analyze: \n {}.".format(lipids_dict)))
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available.".format(lipid))
                # get the metric names available
                metric_names = list(metric[resi][lipid].keys())
                inc = 1
                while len(metric_names) == 0:
                    metric_names = list(metric[resi+inc][lipid].keys())
                    inc += 1
                metric_names_dict = {i: metric_name for i, metric_name in enumerate(metric_names)}
                metric_name = int(input("Please enter the metric name that you want to analyze: \n {}.".format(metric_names_dict)))
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