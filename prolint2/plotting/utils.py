import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib as mpl
from collections import Counter


# define a function to generate a script template
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


# Defining tail strings for different Plotter classes
def point_distribution_tail(name):
    return """
    mean_instance = MeanMetric()
    metric_instance = Metric(contacts, mean_instance)
    mean_contacts = metric_instance.compute()

    # Generate the plot
    PLOT = {}(u, mean_contacts, lipid='CHOL', metric_name='MeanMetric')
    PLOT.save_plot()
            """.format(
        name
    )


def radar_tail(name):
    return """
    metric_instances_list = [MeanMetric(), SumMetric(), MaxMetric()]
    metric_instance = Metric(contacts, metric_instances_list) 
    contacts_out = metric_instance.compute()

    # Generate the plot
    PLOT = {}(contacts_out, resIDs=[2, 3, 5], lipid='POPS', metric_names=['MeanMetric', 'SumMetric', 'MaxMetric'])
    PLOT.save_plot()
            """.format(
        name
    )


def density_map_tail(name):
    return """
    # Generate the plot
    PLOT = {}(u, lipid='CHOL')
    PLOT.save_plot()
            """.format(
        name
    )


def get_metric_list_by_residues(universe, metric, lipid=None, metric_name=None):
    """
    Get a list of metric values for a given lipid and metric name.

    Parameters:
    - universe: Molecular dynamics universe object containing information about the system.
    - metric: Dictionary containing metric data for different residues, lipids, and metrics.
    - lipid: Name of the lipid for which metrics are being extracted (optional).
    - metric_name: Name of the specific metric being extracted (optional).

    Returns:
    - A NumPy array containing the extracted metric values.
    """
    # Extract residue IDs from the universe
    resids = universe.query.residues.resids
    metric_list = []  # List to store extracted metric values

    # Loop through each residue ID
    for resi in resids:
        # Check if the residue ID is present in the metric dictionary
        if resi in metric.keys():
            if metric_name:
                if lipid:
                    try:
                        # Attempt to extract the metric value for the given lipid and metric name
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        # If the metric is not available, append 0 to the list
                        metric_list.append(0)
                else:
                    # Provide user options to select a lipid for analysis
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
                        print("The lipid {} is not available.".format(lipid))
                        break
                else:
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

                # Provide user options to select a metric for analysis
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
            metric_list.append(
                0
            )  # If residue ID is not in the metric dictionary, append 0

    # Convert the list to a NumPy array and return
    return np.array(metric_list)


def get_metrics_for_radar(metrics, metric_names, resIDs=[], lipid=None):
    # Initialize an empty dictionary to store radar metrics for each residue.
    metrics_radar = {}

    # Iterate through the list of residue IDs (resIDs).
    for resi in resIDs:
        # Initialize an empty list for the radar metrics of the current residue.
        metrics_radar[resi] = []

        # Check if the current residue ID (resi) is present in the metrics dictionary.
        if resi in metrics.keys():
            # Iterate through the list of metric names.
            for metric in metric_names:
                # Check if the current metric is present in the metrics dictionary for the given lipid.
                if metric in metrics[resi][lipid].keys():
                    # If the metric is available, add its value to the radar metrics list for the current residue.
                    metrics_radar[resi].append(metrics[resi][lipid][metric])
                else:
                    # If the metric is not available, add 0 to the radar metrics list for the current residue.
                    metrics_radar[resi].append(0)
        else:
            # If the current residue ID is not present in the metrics dictionary,
            # add 0 for all metrics in the radar metrics list for the current residue.
            for metric in metric_names:
                metrics_radar[resi].append(0)

    # Return the radar metrics dictionary and the list of metric names.
    return metrics_radar, metric_names


def compute_density(universe, lipid, size_in_mb):
    # Get the total number of frames in the trajectory
    frames = universe.trajectory.n_frames

    # Create a selection string for the specified lipid
    selection_string = "resname {}".format(lipid)

    # Select the atoms corresponding to the specified lipid
    lipids_ag = universe.select_atoms(selection_string)

    # Loop through the trajectory frames
    for ts in universe.trajectory:
        if ts.frame == 0:
            # Store the positions of the lipid atoms for the first frame
            lipids_xyz = lipids_ag.positions
        else:
            # Append the positions of the lipid atoms for subsequent frames
            lipids_xyz = np.append(lipids_xyz, lipids_ag.positions)

    # Convert the desired size from MB to Bytes
    size = size_in_mb * (1024**2)  # 1 MB = 1024 * 1024 bytes

    # Define a function to slice an array based on a step size
    def slice_array(arr, slice_by):
        mask = np.ones_like(arr, dtype=bool)
        mask[:: int(slice_by)] = False
        xshape = int(arr.shape[0] - mask[:: int(slice_by)].shape[0])
        arr = arr[mask].reshape(xshape, arr.shape[1], 3)
        return arr

    # Reshape the array to represent positions for each frame and atom
    lipids_xyz = lipids_xyz.reshape((frames, lipids_ag.atoms.n_atoms, 3))

    # Reduce the array size to fit within the specified limit
    while lipids_xyz.nbytes > size:
        lipids_xyz = slice_array(lipids_xyz, 10)

    # Reshape the array to a flattened shape for final output
    lipids_xyz = lipids_xyz.reshape(lipids_xyz.shape[0] * lipids_xyz.shape[1], 3)

    # Return the modified array with reduced density of positions
    return lipids_xyz


def shift_range(values, new_min=0, new_max=1):
    """
    Shifts a list of values from their original range to a new range.
    """
    # Find the minimum and maximum values in the input list
    old_min = min(values)
    old_max = max(values)

    # Calculate the range of the original and new ranges
    old_range = old_max - old_min
    new_range = new_max - new_min

    # Initialize a list to store the shifted values
    new_val = []

    # Iterate through the input values and perform the range shifting
    for val in values:
        try:
            # Calculate the new value in the target range using linear mapping
            new_value = (((val - old_min) * new_range) / old_range) + new_min
        except ZeroDivisionError:
            # Handle the case where the range of input values is zero
            new_value = new_min

        # Add the shifted value to the new_val list
        new_val.append(new_value)

    # Return the list of values shifted to the new range
    return new_val


def get_lipid_contact_durations(u, contacts, lipid_type, frequency_filter=0.0):
    """Get the lipid contact frequencies."""

    # Initialize a dictionary to store lipid contact counts
    lipid_contacts = {}

    # Convert residue IDs to a list for easier indexing
    reslist = u.residues.resids.tolist()

    # Loop through each residue in the contact frames
    for resid in contacts.contact_frames:
        for lipid in contacts.contact_frames[resid]:
            # Check if the lipid's resname matches the desired lipid_type
            if u.residues.resnames[reslist.index(lipid)] == lipid_type:
                # If lipid not encountered before, add it to the dictionary
                if lipid not in lipid_contacts.keys():
                    lipid_contacts[lipid] = 0

                # Increment the contact count for the lipid
                lipid_contacts[lipid] += len(contacts.contact_frames[resid][lipid])

    # Calculate the total number of frames (trajectory snapshots)
    n_frames = u.trajectory.n_frames - 1

    # Filter lipid_contacts based on the specified frequency_filter
    lipid_contacts = {
        lipid: frames / n_frames
        for lipid, frames in lipid_contacts.items()
        if frames / n_frames > frequency_filter
    }

    # Sort the lipid_contacts by frequency in descending order
    temp = sorted(lipid_contacts.items(), key=lambda x: x[1], reverse=True)

    # Create a DataFrame to store the sorted lipid IDs and their frequencies
    import pandas as pd

    result_df = pd.DataFrame(
        {"Lipip ID": [i[0] for i in temp], "Frequency": [i[1] for i in temp]}
    )

    # Return the DataFrame containing lipid contact information
    return result_df


def inverse_dict_keys(d):
    """
    Takes a dictionary as input and returns a new dictionary where the keys of the input dictionary
    are reversed in order while preserving their corresponding values.
    """
    reversed_keys = list(d.keys())[
        ::-1
    ]  # Get the keys of the input dictionary in reversed order
    inverted_dict = {
        key: d[key] for key in reversed_keys
    }  # Create a new dictionary with reversed keys and original values
    return inverted_dict  # Return the resulting inverted dictionary


def create_logo_df(universe, metric, **kwargs):
    # Extract residue IDs from the provided universe
    resids = universe.query.residues.resids

    # Convert amino acid codes to full names for each residue
    resnames = [
        mda.lib.util.convert_aa_code(x) for x in universe.query.residues.resnames
    ]

    # Calculate the specified metric values for each residue using provided kwargs
    res_metrics = get_metric_list_by_residues(universe, metric, **kwargs)

    # Create a pandas DataFrame to store the logo data
    df = pd.DataFrame(
        {
            "ResID": resids,  # Column for residue IDs
            "Resname": resnames,  # Column for residue names
            "Occupancy": [0.75]
            * len(resids),  # Default occupancy value for each residue
            "Metric": res_metrics,  # Column for calculated metric values
        }
    )

    return df  # Return the created DataFrame containing the logo data


def get_frame_contact_intervals(frames, continuity_filter, tolerance):
    """
    Get the intervals of frames in which a contact is present.

    Args:
        frames (list): A list of frames in which a contact is present.
        continuity_filter (int): Minimum duration of a contact interval to be considered valid.
        tolerance (int): The number of frames to tolerate before considering a new interval.

    Returns:
        ranges_collect (list): A list of tuples containing the start and end frames of each
            valid interval.

    """

    # Initialize an empty list to collect the intervals
    ranges_collect = []

    # Initialize a variable to track the start of the current interval
    range_start = 0

    # Loop through each frame and its index in the list
    for ix, el in enumerate(frames):
        if ix == 0:
            range_start = el
            continue

        # Get the previous frame
        prev_el = frames[ix - 1]

        # Check if the current frame is outside the tolerance range from the previous frame
        if not el - tolerance <= prev_el:
            # The previous interval is complete, add it to the collection
            ranges_collect.append((range_start, prev_el))
            # Start a new interval from the current frame
            range_start = el

        # Check if this is the last frame in the list
        if ix == len(frames) - 1:
            # Add the last interval to the collection
            ranges_collect.append((range_start, el))

    # Filter and return the collected intervals based on continuity_filter
    return [
        (pair[0], pair[1])
        for pair in ranges_collect
        if pair[1] - pair[0] >= continuity_filter
    ]


# Define a function to reverse the colors of a given colormap
def reverse_colourmap(cmap, name="my_cmap_r"):
    # Initialize an empty list to store the reversed color data
    reverse = []
    # Initialize an empty list to store the keys of the color channels
    k = []

    # Loop through each key in the color map's segment data
    for key in cmap._segmentdata:
        # Append the key to the list of keys
        k.append(key)
        # Retrieve the color data for the current channel
        channel = cmap._segmentdata[key]
        # Initialize an empty list to store the reversed color data for this channel
        data = []

        # Loop through each color point in the channel
        for t in channel:
            # Reverse the first value (position) and swap the second and third values (green and blue)
            data.append((1 - t[0], t[2], t[1]))

        # Add the reversed color data for this channel to the list
        reverse.append(sorted(data))

    # Create a dictionary that maps keys to their corresponding reversed color data
    LinearL = dict(zip(k, reverse))
    # Create a new LinearSegmentedColormap using the reversed color data
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)

    # Return the newly created reversed colormap
    return my_cmap_r


def get_lipid_contact_frequencies(u, contacts, lipid_type):
    """Get the lipid contact frequencies."""

    # Nested function to sort a contact list by frequency
    def sort_by_frequency(contact_list):
        contact_list.sort(key=lambda x: x[1], reverse=True)
        return contact_list

    # Threshold for excluding contacts with low frequency
    contact_threshold = 0

    # Initialize dictionaries to store lipid and residue contact frequencies
    lipid_frequency = {lipid: {} for lipid in u.database.unique_resnames}
    residue_contact_freq = {}

    # Loop through residues and their lipid contacts
    for residue, lipid_contacts in contacts.compute_metric("sum").items():
        for lipid, contact_counter in lipid_contacts.items():
            # Sort the contact_counter dictionary by its values (frequencies)
            sorted_contacts = sorted(
                contact_counter.items(), key=lambda x: x[1], reverse=True
            )

            # Iterate over sorted contacts
            for lipid_id, freq in sorted_contacts:
                # Exclude contacts with frequency below the threshold
                if freq <= contact_threshold:
                    continue

                # Update lipid_frequency dictionary
                if lipid_id in lipid_frequency[lipid]:
                    lipid_frequency[lipid][int(lipid_id)] += freq
                else:
                    lipid_frequency[lipid][int(lipid_id)] = freq

                # Update residue_contact_freq dictionary
                if int(lipid_id) in residue_contact_freq:
                    residue_contact_freq[int(lipid_id)].append((int(residue), freq))
                else:
                    residue_contact_freq[int(lipid_id)] = [(int(residue), freq)]

    # Sort and store lipid contact frequencies
    for lipid, values in lipid_frequency.items():
        lipid_frequency[lipid] = Counter(values).most_common()

    # Sort and store residue contact frequencies
    for lipid_id, vals in residue_contact_freq.items():
        residue_contact_freq[lipid_id] = sort_by_frequency(vals)

    # Return the sorted lipid contact frequencies and residue contact frequencies for a given lipid type
    return [x[0] for x in lipid_frequency[lipid_type]], residue_contact_freq
