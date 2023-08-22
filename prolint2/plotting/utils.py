import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib as mpl

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
from prolint2.plotting import Plotter
from prolint2.plotting.utils import *

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


def compute_density(universe, lipid, size_in_mb):
    frames = universe.trajectory.n_frames
    selection_string = "resname {}".format(lipid)

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


def shift_range(values, new_min=0, new_max=1):
    """
    Shift values from one range to another.
    """
    old_min = min(values)
    old_max = max(values)

    old_range = old_max - old_min
    new_range = new_max - new_min
    new_val = []
    for val in values:
        try:
            new_value = (((val - old_min) * new_range) / old_range) + new_min
        except ZeroDivisionError:
            new_value = new_min
        new_val.append(new_value)
    return new_val


def get_lipid_contact_frequencies(u, contacts, lipid_type, frequency_filter=0.0):
    """Get the lipid contact frequencies."""
    lipid_contacts = {}
    reslist = u.residues.resids.tolist()

    for resid in contacts.contact_frames:
        for lipid in contacts.contact_frames[resid]:
            # filter by resname
            if u.residues.resnames[reslist.index(lipid)] == lipid_type:
                if lipid not in lipid_contacts.keys():
                    lipid_contacts[lipid] = 0
                lipid_contacts[lipid] += len(contacts.contact_frames[resid][lipid])

    n_frames = u.trajectory.n_frames - 1

    lipid_contacts = {
        lipid: frames / n_frames
        for lipid, frames in lipid_contacts.items()
        if frames / n_frames > frequency_filter
    }
    temp = sorted(lipid_contacts.items(), key=lambda x: x[1], reverse=True)
    return pd.DataFrame(
        {"Lipip ID": [i[0] for i in temp], "Frequency": [i[1] for i in temp]}
    )


def inverse_dict_keys(d):
    reversed_keys = list(d.keys())[::-1]
    inverted_dict = {key: d[key] for key in reversed_keys}
    return inverted_dict


def create_logo_df(universe, metric, **kwargs):
    # create pandas dataframe
    resids = universe.query.residues.resids
    resnames = [
        mda.lib.util.convert_aa_code(x) for x in universe.query.residues.resnames
    ]
    res_metrics = get_metric_list_by_residues(universe, metric, **kwargs)
    df = pd.DataFrame(
        {
            "ResID": resids,
            "Resname": resnames,
            "Occupancy": [0.75] * len(resids),
            "Metric": res_metrics,
        }
    )
    return df


def get_frame_contact_intervals(frames, continuity_filter, tolerance):
    """
    Get the intervals of frames in which a contact is present.

    Args:
        frames (list): A list of frames in which a contact is present.
        tolerance (int): The number of frames to tolerate before considering a new interval.

    Returns:
        ranges_collect (list): A list of tuples containing the start and end frames of each
            interval.

    """
    ranges_collect = []
    range_start = 0
    for ix, el in enumerate(frames):
        if ix == 0:
            range_start = el
            continue

        prev_el = frames[ix - 1]
        if not el - tolerance <= prev_el:
            ranges_collect.append((range_start, prev_el))
            range_start = el
        if ix == len(frames) - 1:
            ranges_collect.append((range_start, el))
    return [
        (pair[0], pair[1])
        for pair in ranges_collect
        if pair[1] - pair[0] >= continuity_filter
    ]


def reverse_colourmap(cmap, name = 'my_cmap_r'):    
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r