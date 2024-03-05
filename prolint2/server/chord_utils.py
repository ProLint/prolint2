#!/usr/bin/env python

import numpy as np
from itertools import combinations
from prolint2.server.utils import calculate_contact_intervals

import os
import configparser

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini"))
parameters_config = config["Parameters"]

# def per_lipid_contacts(ts, lipids, frame_cutoff=10):
#     """
#     Given a list of lipid IDs, returns a dict with these lipid IDs
#     as keys, and values set to a dict containing the residues these lipids
#     interact with as keys and the corresponding number of contacts as values.
#     These contacts can be filtered using the `frame_cutoff` option.

#     TODO:
#     `frame_cutoff` should operate on a percentage of trajectory length.
#     """
#     results = {k: {} for k in lipids}
#     for k, v in ts.contacts.contact_frames.items():
#         if len(v) < frame_cutoff:
#             continue
#         # r, l = [int(x) for x in k.split(',')] # k used to be a string formatted as 'residue,lipid'
#         # k now is a tuple of (residue, lipid)
#         r, l = k
#         if l in lipids:
#             results[l][r] = len(v)
#     return results


def per_lipid_contacts(ts, lipids, frame_cutoff=10):
    results = {k: {} for k in lipids}
    for residue_id, lipid_dict in ts.contacts.contact_frames.items():
        for lipid_id, frames in lipid_dict.items():
            if len(frames) < frame_cutoff:
                continue
            if lipid_id in lipids:
                results[lipid_id][residue_id] = len(frames)
    return results


def sort_dict(d, cutoff=None):
    """
    Takes a dictionary as input, and sorts it according to values.
    The purpose of this function is to order residue:contacts dicts build
    by `per_lipid_contacts` function according to the number of contacts.

    The `cutoff` value returns only a subsection of this dict. The cutoff
    is used by `get_linked_nodes` to get a subsection of edges/links. This is
    important because visually we can not easily distinguish between link widths,
    so filtering here allows us to only show the top X links.
    """
    item_list = sorted(d.items(), key=lambda x: x[1], reverse=True)
    if cutoff:
        return dict(item_list[:cutoff])
    return dict(item_list)


def residue_pair_matching_contacts(res1_contacts, res2_contacts):
    """
    docstring: important improvement to the network app contacts
    """
    total_contacts = 0
    for s1, e1 in res1_contacts:
        for s2, e2 in res2_contacts:
            if not e1 < s2 and not s1 > e2:
                contacts = min(e1, e2) - max(s1, s2)
                total_contacts += contacts

    return total_contacts


def get_ordered_combinations(lipid_contacts):
    """
    This function turns residue:contact information into shared contacts. To speed up
    performance, shared contacts are define in an adhoc manner. It takes as input the
    returned dict by `per_lipid_contacts`. It then forms all possible size-2 combinations.

    For a lipid ID, all these residue-residue combinations are seen as shared contacts.
    That is, if a lipid interacts with residues 1, 2, and 3, then all possible size-2
    combination (1-2, 1-3, 2-3) are seen as shared contacts.

    These combinations are assigned as keys, with values corresponding to the sum
    of all matching shared contacts across all lipid IDs. This dict is then returned.
    This approach seems to work quite well in addition to being intuitive. One potential
    disadvantage, however, is that if a residue has very high contact count, then so will
    all combinations involving that residue.

    Another important point to consider is that since shared contacts are defined adhoc,
    an edge linking two nodes, cannot reliably be interpreted as geometric distance.
    Meaning, it does not mean that those two nodes/residues are connected by interactions
    with the same lipid.

    TODO:
    1.
    Is it possible to improve the last point, by counting only real shared contacts?
    The information is already contained in `contacts.contact_frames`.
    2.
    combinations call already returns an iterable. Simply by removing the
    subsequent call to list().
    """
    ordered_combinations = {}
    for lipid_id, lipid_vals in lipid_contacts.items():
        sorted_lipid_vals = sort_dict(lipid_vals)

        lipid_residues = list(sorted_lipid_vals.keys())
        for res1, res2 in list(combinations(lipid_residues, 2)):
            value = sum([sorted_lipid_vals[res1], sorted_lipid_vals[res2]])
            key = f"{res1},{res2}"

            if key in ordered_combinations:
                ordered_combinations[key] = ordered_combinations[key] + value
            else:
                ordered_combinations[key] = value

    return ordered_combinations


def shared_contacts(contacts, top_lipids, lipid_contact_frames, *args, **kwargs):
    """
    Aim: improve the shortcomings outlined in `get_ordered_combinations`.
    """
    lipid_shared_contacts = {}
    for lipid in top_lipids:
        contact_intervals = calculate_contact_intervals(
            contacts, lipid_contact_frames, lipid, *args, **kwargs
        )
        residue_contacts = {}
        for res1, res2 in combinations(contact_intervals.keys(), 2):
            pair_contacts = residue_pair_matching_contacts(
                contact_intervals[res1], contact_intervals[res2]
            )
            residue_contacts[f"{res1},{res2}"] = pair_contacts
        lipid_shared_contacts[lipid] = residue_contacts

    shared_contacts_all = {}
    for residue_pair_contacts in lipid_shared_contacts.values():
        for residue_pair, pair_contacts in residue_pair_contacts.items():
            if residue_pair in shared_contacts_all:
                shared_contacts_all[residue_pair] += pair_contacts
            else:
                shared_contacts_all[residue_pair] = pair_contacts

    # return sort_dict(shared_contacts_all)
    return lipid_shared_contacts, sort_dict(shared_contacts_all)


def get_linked_nodes(ordered_combinations, cutoff=100):
    """
    This function along with `get_node_list` help deal with amChart shortcomings
    related to ChordNonRibbon. We need to add hidden nodes to keep the chord
    architecture nice and matching the structure of the protein. Here we get
    the node ids that have data linked to them.
    """
    linked_nodes = []
    for residue_key in sort_dict(ordered_combinations, cutoff=cutoff).keys():
        res1, res2 = [int(res) for res in residue_key.split(",")]
        linked_nodes.extend([res1, res2])
    return np.unique(linked_nodes).tolist()


def get_node_list(n_residues, linked_nodes):
    """
    This function along with `get_linked_nodes` help deal with amChart shortcomings
    related to ChordNonRibbon. We need to add hidden nodes to keep the chord
    architecture nice and matching the structure of the protein. Here we get
    the node ids that do NOT have data linked to them.

    TODO:
    We need to have checks on the total number of nodes that can be visualized.
    Currently, we add 200 nodes on top of the already existing linked nodes.
    """

    hidden_nodes = 200

    node_list = np.linspace(0, n_residues, hidden_nodes, dtype=int)
    node_list = [int(x) for x in node_list if int(x) not in linked_nodes]
    nodes = sorted(node_list + linked_nodes)

    hidden_node_indices = [ix for (ix, x) in enumerate(nodes) if x not in linked_nodes]
    return nodes, hidden_node_indices


def get_chord_elements(ts, nodes, ordered_combinations, cutoff=500):
    """
    Prepares the input data so it can be read and understood by the amCharts.
    The output is a simple list of dicts. Each entry containing info for
    one link. We also add residue ID and residue name information which is
    important to link node data with other apps.
    """
    node_links = list(combinations(nodes, 2))
    position_node_links = [x for x in node_links if x[0] == 0]

    resnums = ts.query.residues.resnums
    resnames = ts.query.residues.resnames
    node_names = {x[0]: f"{x[0]} {x[1]}" for x in list(zip(resnums, resnames))}

    chord_elements = []
    for res1, res2 in position_node_links:
        chord_elements.append({"from": 0, "to": node_names[res2], "value": 0})

    contact_max = max(ordered_combinations.values())
    for k, v in sort_dict(ordered_combinations, cutoff=cutoff).items():
        res1, res2 = k.split(",")
        # This will break the nodes. We need to add the nodes we are skipping here
        # to the position_residue_indices
        # if abs(int(res1) - int(res2)) < 4:
        #     continue

        chord_elements.append(
            {
                "from": node_names[int(res1)],
                "to": node_names[int(res2)],
                "value": v,
                "valueWidth": float(v) / contact_max,
            }
        )

    return chord_elements


def contact_chord(ts, contacts, top_lipid_ids, lipid_contact_frames, cutoff=100):
    """
    We call all functions here. We return the chord elements (these are the data
    amCharts needs to render nodes and links), we also return information on which
    nodes to hide, and which nodes correspond to which lipid ID.

    TODO:
    `per_lipid_nodes` can be further refined.
    """
    # lipid_contacts = per_lipid_contacts(ts, top_lipid_ids)
    lipid_shared_contacts, ordered_combinations = shared_contacts(
        contacts,
        top_lipid_ids,
        lipid_contact_frames,
        residues_to_show=int(parameters_config["residues_to_show"]),
        intervals_to_filter_out=int(parameters_config["intervals_to_filter_out"]),
    )
    # ordered_combinations = get_ordered_combinations(lipid_contacts)
    linked_nodes = get_linked_nodes(ordered_combinations, cutoff=cutoff)
    nodes, hidden_node_indices = get_node_list(ts.query.n_residues, linked_nodes)
    chord_elements = get_chord_elements(ts, nodes, ordered_combinations, cutoff=cutoff)

    per_lipid_nodes = {}
    for lipid, pair_contacts in lipid_shared_contacts.items():
        lipid_nodes = get_linked_nodes(pair_contacts, cutoff=None)
        all_lipid_nodes = [x for x in lipid_nodes if x in linked_nodes]
        if all_lipid_nodes:
            per_lipid_nodes[lipid] = all_lipid_nodes

    return chord_elements, hidden_node_indices, per_lipid_nodes
