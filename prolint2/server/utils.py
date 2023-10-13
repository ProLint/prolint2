#!/usr/bin/env python

import os
import configparser

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini"))
parameters_config = config["Parameters"]


def get_frame_contact_intervals(frames, tolerance=int(parameters_config["tolerance"])):
    """
    Get frame ranges.

    TODO:
    write doc
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
    return ranges_collect


def calculate_contact_intervals(
    contacts,
    g,
    lipid_id,
    residues_to_show=int(parameters_config["residues_to_show"]),
    intervals_to_filter_out=int(parameters_config["intervals_to_filter_out"]),
):
    """
    TODO:
    write doc
    """
    contact_intervals = {}
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = contacts.contact_frames[res][lipid_id]
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue

            if res in contact_intervals:
                contact_intervals[res].append((start, end))
            else:
                contact_intervals[res] = [(start, end)]

    return contact_intervals


def amCharts_contact_intervals(
    contacts,
    g,
    lipid_id,
    residues_to_show=int(parameters_config["residues_to_show"]),
    intervals_to_filter_out=int(parameters_config["intervals_to_filter_out"]),
):
    """
    TODO:
    write doc
    """
    contact_intervals = []
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = contacts.contact_frames[res][lipid_id]
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue

            d = {
                "category": res,
                "startFrame": start,
                "endFrame": end,
                "lipid_id": lipid_id,
            }
            contact_intervals.append(d)

    return contact_intervals
