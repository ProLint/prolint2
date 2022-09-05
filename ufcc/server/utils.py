#!/usr/bin/env python

def get_frame_contact_intervals(frames, tolerance=6):
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

        prev_el = frames[ix-1]
        if not el-tolerance <= prev_el:
            ranges_collect.append((range_start, prev_el))
            range_start = el
        if ix == len(frames) - 1:
            ranges_collect.append((range_start, el))
    return ranges_collect

def calculate_contact_intervals(TS, g, lipid_id, residues_to_show=15, intervals_to_filter_out=10):
    """
    TODO:
    write doc
    """
    contact_intervals = {}
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = TS.contacts.contact_frames[f'{res},{lipid_id}']
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue

            if res in contact_intervals:
                contact_intervals[res].append((start, end))
            else:
                contact_intervals[res] = [(start, end)]

    return contact_intervals

def amCharts_contact_intervals(TS, g, lipid_id, residues_to_show=15, intervals_to_filter_out=10):
    """
    TODO:
    write doc
    """
    contact_intervals = []
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = TS.contacts.contact_frames[f'{res},{lipid_id}']
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue

            d = {
            "category": res,
            "startFrame": start,
            "endFrame": end,
            "lipid_id": lipid_id
            }
            contact_intervals.append(d)

    return contact_intervals
