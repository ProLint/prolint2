#!/usr/bin/env python

import numpy as np
from itertools import combinations

def per_lipid_contacts(ts, lipids, frame_cutoff=10):
    results = {k: {} for k in lipids}
    for k, v in ts.contacts.contact_frames.items():
        if len(v) < frame_cutoff:
            continue
        r, l = [int(x) for x in k.split(',')]
        if l in lipids:
            results[l][r] = len(v)
    return results

def sort_dict(d, cutoff=None):
    item_list = sorted(d.items(), key=lambda x: x[1], reverse=True)
    if cutoff:
        return dict(item_list[:cutoff])
    return dict(item_list)

def get_ordered_combinations(lipid_contacts):
    ordered_combinations = {}
    for lipid_id, lipid_vals in lipid_contacts.items():
        sorted_lipid_vals = sort_dict(lipid_vals)

        lipid_residues = list(sorted_lipid_vals.keys())
        for res1, res2 in list(combinations(lipid_residues, 2)):
            value = sum([sorted_lipid_vals[res1], sorted_lipid_vals[res2]])
            key = f'{res1},{res2}'

            if key in ordered_combinations:
                ordered_combinations[key] = ordered_combinations[key] + value
            else:
                ordered_combinations[key] = value

    return ordered_combinations

def get_linked_nodes(ordered_combinations, cutoff=100):
    linked_nodes = []
    for residue_key in sort_dict(ordered_combinations, cutoff=cutoff).keys():
        res1, res2 = [int(res) for res in residue_key.split(',')]
        linked_nodes.extend([res1, res2])
    return np.unique(linked_nodes).tolist()

def get_node_list(n_residues, linked_nodes):
    hidden_nodes = 200

    node_list = np.linspace(0, n_residues, hidden_nodes, dtype=int)
    node_list = [int(x) for x in node_list if int(x) not in linked_nodes]
    nodes = sorted(node_list + linked_nodes)

    hidden_node_indices = [ix for (ix, x) in enumerate(nodes) if x not in linked_nodes]
    return nodes, hidden_node_indices

def get_chord_elements(ts, nodes, ordered_combinations, cutoff=500):
    node_links = list(combinations(nodes, 2))
    position_node_links = [x for x in node_links if x[0] == 0]

    resnums = ts.query.selected.residues.resnums
    resnames = ts.query.selected.residues.resnames
    node_names = {x[0]: f'{x[0]} {x[1]}' for x in list(zip(resnums, resnames))}

    chord_elements = []
    for res1, res2 in position_node_links:
        chord_elements.append({
            "from": 0,
            "to": node_names[res2],
            "value": 0
        })

    contact_max = max(ordered_combinations.values())
    for k, v in sort_dict(ordered_combinations, cutoff=cutoff).items():
        res1, res2 = k.split(',')
        # This will break the nodes. We need to add the nodes we are skipping here
        # to the position_residue_indices
        # if abs(int(res1) - int(res2)) < 4:
        #     continue

        chord_elements.append({
            "from": node_names[int(res1)],
            "to": node_names[int(res2)],
            "value": v,
            "valueWidth": float(v) / contact_max
        })

    return chord_elements

def contact_chord(ts, top_lipid_ids, cutoff=100):
    lipid_contacts = per_lipid_contacts(ts, top_lipid_ids)
    ordered_combinations = get_ordered_combinations(lipid_contacts)
    linked_nodes = get_linked_nodes(ordered_combinations, cutoff=cutoff)
    nodes, hidden_node_indices = get_node_list(ts.query.selected.n_residues, linked_nodes)
    chord_elements = get_chord_elements(ts, nodes, ordered_combinations, cutoff=cutoff)

    return chord_elements, hidden_node_indices