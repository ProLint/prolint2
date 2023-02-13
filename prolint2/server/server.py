from collections import Counter
from prolint2.interactive_sel import interactive_selection
import os
import ast
import json

from bottle import route, run, template, debug, static_file, request
from prolint2.contacts import SerialDistances

import MDAnalysis as mda
from prolint2.prolint2 import PL2
from io import StringIO

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))

BACKEND_DATA = None
TS = None
ARGS = None
data = None
data_loaded = False


def sort_lipids(ts):
    """
    Sort lipid contacts according to their contact frequency, all the while keeping track
    of residue IDs, and number of contacts with each residue.

    Returns:
        t (dict): Stores lipid IDs and their contact frequency, sorted in descending order
        g (dict): For each lipid ID, stores the residues in contact and the corresponding
                  frequency.
    """

    def sort_tuple(tup):
        tup.sort(key=lambda x: x[1], reverse=True)
        return tup

    # TODO:
    # top lipid number should be put in the config.
    contact_threshold = ts.n_frames * 0.05

    # initialize dictionary to store values:
    t = {k: {} for k in ts.database_unique}
    g = {}
    for ix, (residue, lipid_contacts) in enumerate(ts.contacts.contacts.items()):
        for lipid, contact_counter in lipid_contacts.items():
            top10_counter = contact_counter.most_common()
            for (lipid_id, lipid_counter) in top10_counter:
                # Exclude short-lived contacts
                if lipid_counter <= contact_threshold:
                    continue
                if lipid_id in t[lipid]:
                    t[lipid][lipid_id] += lipid_counter
                    g[lipid_id].append((residue, lipid_counter))
                else:
                    t[lipid][lipid_id] = lipid_counter
                    g[lipid_id] = [(residue, lipid_counter)]

    for lipid, values in t.items():
        t[lipid] = Counter(values).most_common()

    # for lipid, values in g.items():
    for lipid_id, vals in g.items():
        g[lipid_id] = sort_tuple(vals)

    return t, g


def get_frame_contact_intervals(frames, tolerance=6):
    """
    Get frame ranges
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


def get_gantt_app_data(g, lipid_id, residues_to_show=15, intervals_to_filter_out=10):
    gantt_data = []
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = TS.contacts.contact_frames[f"{res},{lipid_id}"]
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue
            gantt_data.append(
                {
                    # "category": f'{res}',
                    "category": res,
                    "startFrame": start,
                    "endFrame": end,
                    "lipid_id": lipid_id,
                }
            )

    categories = []
    for y in [x["category"] for x in gantt_data]:
        if y not in categories:
            categories.append(y)
    return gantt_data, categories


@route("/static/<filepath:path>")
def server_static(filepath):
    return static_file(filepath, root=os.path.join(SERVER_PATH, "static"))


@route("/")
def index():
    return template(os.path.join(SERVER_PATH, "home.tpl"))


@route("/app")
def app():
    return static_file("index.html", root=SERVER_PATH)


@route("/prolint2")
def prolint2():
    import sys

    print(request.body.getvalue().decode("utf-8"), file=sys.stdout)
    return request.body


@route("/toplipids/:metadata")
def top_lipid_listener(metadata):
    global BACKEND_DATA

    metadata = ast.literal_eval(metadata)
    lipid_id = metadata["lipidID"]

    gantt_data, categories = get_gantt_app_data(
        BACKEND_DATA["lipid_contact_frames"], lipid_id
    )
    # This will sort the residues
    # sorted_gantt_data = sorted(gantt_data, key=lambda d: d['category'])

    # ags = TS.query.selected.select_atoms(f'resid {" ".join([str(x) for x in categories])}')
    # labeled_categories = [[int(x.resid), x.resname] for x in ags]
    return {
        "ganttData": gantt_data,
        "topLipids": categories,
    }


@route("/distance/:metadata")
def distance_array_listener(metadata):
    global BACKEND_DATA
    global TS

    metadata = ast.literal_eval(metadata)
    lipid_id = metadata["lipidID"]
    residue_id = int(metadata["residueID"])

    ri = SerialDistances(
        TS.query.selected.universe,
        TS.query.selected,
        TS.database.selected,
        lipid_id,
        residue_id,
        TS.contacts.contact_frames[f"{residue_id},{lipid_id}"],
    )
    ri.run(verbose=False)

    hm_data, la_data = [], []
    for lx, la in enumerate(ri.lipid_atomnames):
        la_data.append({"LipidAtoms": la})
        for rx, ra in enumerate(ri.resid_atomnames):
            v = ri.distance_array[lx, rx]
            hm_data.append({"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)})
    ra_data = [{"ResidueAtoms": x} for x in ri.resid_atomnames]

    return {
        "heatmapData": hm_data,
        "lipidAtomsData": la_data,
        "residueAtomsData": ra_data,
    }


@route("/tabledata/:metadata")
def table_listener(metadata):
    global BACKEND_DATA

    metadata = ast.literal_eval(metadata)
    lipid = metadata["lipid"]

    table_data = []
    for ix, (lipid_id, freq) in enumerate(BACKEND_DATA["top_lipids"][lipid]):
        table_data.append({"id": ix, "lipidID": lipid_id, "contactFrequency": freq})

    return {
        "tableData": table_data,
    }


@route("/pdb/:metadata")
def blob(metadata):

    global ARGS
    u = mda.Universe(ARGS.structure, ARGS.trajectory)
    protein = u.select_atoms("protein")
    pstream = mda.lib.util.NamedStream(StringIO(), "dummy.pdb")
    with mda.Writer(pstream, format="PDB") as w:
        w.write(protein)

    return pstream.read()


@route("/data/:metadata")
def listener(metadata):

    global data_loaded
    global data
    global BACKEND_DATA
    global TS

    # TODO:
    # Bottle should provide the metadata already,
    # perhaps via the following:
    # from bottle import response, request
    metadata = ast.literal_eval(metadata)

    lipid = metadata["lipid"]
    protein = metadata["protein"]

    if lipid == "" and protein == "":
        # Starting setup:
        try:
            lipid = BACKEND_DATA["lipids"][0]
            protein = BACKEND_DATA["proteins"][0]
        except:
            print("Detached EXECUTION")
            print("This is currently meant for testing only. Not guaranteed to work!")
            BACKEND_DATA = independent_execution()
            lipid = BACKEND_DATA["lipids"][0]
            protein = BACKEND_DATA["proteins"][0]

    table_data = []
    for ix, (lipid_id, freq) in enumerate(BACKEND_DATA["top_lipids"][lipid]):
        table_data.append({"id": ix, "lipidID": lipid_id, "contactFrequency": freq})

    # Initiate ganttApp with the top lipid
    lipid_id = BACKEND_DATA["top_lipids"][lipid][0][0]
    gantt_data, categories = get_gantt_app_data(
        BACKEND_DATA["lipid_contact_frames"], lipid_id
    )

    # Initiate heatmapApp with the top residue
    residue_id = BACKEND_DATA["lipid_contact_frames"][lipid_id][0][0]
    ri = SerialDistances(
        TS.query.selected.universe,
        TS.query.selected,
        TS.database.selected,
        lipid_id,
        residue_id,
        TS.contacts.contact_frames[f"{residue_id},{lipid_id}"],
    )
    ri.run(verbose=False)

    hm_data, la_data = [], []
    for lx, la in enumerate(ri.lipid_atomnames):
        la_data.append({"LipidAtoms": la})
        for rx, ra in enumerate(ri.resid_atomnames):
            v = ri.distance_array[lx, rx]
            hm_data.append({"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)})
    ra_data = [{"ResidueAtoms": x} for x in ri.resid_atomnames]

    # TODO:
    # Possibly, avoid single point of failure on these dictionary lookups?
    response = {
        "data": BACKEND_DATA["data"][protein][lipid],
        "proteins": BACKEND_DATA["proteins"],
        "lipids": BACKEND_DATA["lipids"],
        "pieData": BACKEND_DATA["pie_data"],
        "ganttData": gantt_data,
        "topLipids": categories,
        "globalTopLipids": BACKEND_DATA["top_lipids"],
        "lipidContactFrames": BACKEND_DATA["lipid_contact_frames"],
        "tableData": table_data,
        "heatmapData": hm_data,
        "lipidAtomsData": la_data,
        "residueAtomsData": ra_data,
        "frameNumber": TS.n_frames,
    }
    return response


def start_server(
    payload=None, debug_bool=False, reloader=True, port=8351, i_bool=True, e_file=False
):

    global ARGS
    # ProLint2 calls:
    args = payload
    ARGS = args
    ts = PL2(args.structure, args.trajectory, add_lipid_types=args.other_lipids)
    # For interactive selection of the groups for the contacts calculation
    if i_bool:
        ts = interactive_selection(ts)
    ts.contacts.compute(cutoff=args.cutoff)

    # for exporting the data
    if e_file:
        ts.contacts.export(args.e_file)

    payload = ts.contacts.server_payload()

    t, g = sort_lipids(ts)
    payload["top_lipids"] = t
    payload["lipid_contact_frames"] = g

    # Make data accessible globally
    global BACKEND_DATA
    global TS
    BACKEND_DATA = payload
    TS = ts

    debug(debug_bool)
    run(reloader=reloader, host="localhost", port=port)


def independent_execution():
    """
    If we are not calling the server through the prolint2 executable, but
    independently, locally for testing purposes, we will load local data file
    and serve that to the dashboard.
    """
    with open(os.path.join(SERVER_PATH, "girk.json"), "r") as fp:
        data = json.load(fp)

    pie_data = [
        {
            "category": "LocalGirk",
            "value": 500,
            "subData": [
                {"category": "CHOL", "value": 300},
                {"category": "POPE", "value": 150},
                {"category": "POPS", "value": 50},
            ],
        }
    ]

    # ganttApp data input requirement
    gantt_data = [
        {
            "category": "Lipid 1",
            "startFrame": 0,
            "endFrame": 10,
        },
        {
            "category": "Lipid 1",
            "startFrame": 45,
            "endFrame": 75,
        },
        {
            "category": "Lipid 1",
            "startFrame": 90,
            "endFrame": 100,
        },
        {
            "category": "Lipid 2",
            "startFrame": 10,
            "endFrame": 35,
        },
        {
            "category": "Lipid 2",
            "startFrame": 45,
            "endFrame": 60,
        },
    ]
    top_10_lipids = ["Lipid 1", "Lipid 2"]

    payload = {
        "data": data,
        "proteins": ["LocalGirk"],
        "lipids": list(data["LocalGirk"].keys()),
        "pie_data": pie_data,
        "gantt_data": gantt_data,
        "top_10_lipids": top_10_lipids,
    }

    return payload


if __name__ == "__main__":
    start_server(debug_bool=True)
