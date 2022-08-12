from collections import Counter
import os
import ast
import json

from bottle import route, run, template, debug, static_file, request
from ufcc.contacts import ProLintSerialDistances

from ufcc.ufcc import UFCC

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))

BACKEND_DATA = None
TS = None
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
        tup.sort(key = lambda x: x[1], reverse=True)
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
                if lipid_counter <= contact_threshold: continue
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

        prev_el = frames[ix-1]
        if not el-tolerance <= prev_el:
            ranges_collect.append((range_start, prev_el))
            range_start = el
        if ix == len(frames) - 1:
            ranges_collect.append((range_start, el))
    return ranges_collect

def get_gantt_app_data(g, lipid_id, residues_to_show=15, intervals_to_filter_out=10):
    gantt_data = []
    for res, _ in g[lipid_id][:residues_to_show]:
        frame_numbers = TS.contacts.contact_frames[f'{res},{lipid_id}']
        frame_intervals = get_frame_contact_intervals(frame_numbers)
        for start, end in frame_intervals:
            if end - start < intervals_to_filter_out:
                continue
            gantt_data.append({
            "category": f'{res}',
            "startFrame": start,
            "endFrame": end,
            })

    categories = []
    for x in [x['category'] for x in gantt_data]:
        if x not in categories:
            categories.append(x)
    return gantt_data, categories


@route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root=os.path.join(SERVER_PATH, 'static'))

@route('/')
def index():
    return template(os.path.join(SERVER_PATH, 'home.tpl'))

@route('/app')
def app():
    return static_file('index.html', root=SERVER_PATH)

@route('/ufcc')
def ufcc():
    import sys
    print(request.body.getvalue().decode('utf-8'), file=sys.stdout)
    return request.body

@route('/toplipids/:metadata')
def top_lipid_listener(metadata):
    global BACKEND_DATA

    metadata = ast.literal_eval(metadata)
    lipid_id = metadata['lipidID']
    gantt_data, categories = get_gantt_app_data(BACKEND_DATA['lipid_contact_frames'], lipid_id)
    print ('cat length', len(categories))

    return {
        "ganttData": gantt_data,
        "topLipids": categories,
    }

@route('/data/:metadata')
def listener(metadata):

    global data_loaded
    global data
    global BACKEND_DATA
    global TS

    # print ('TS object: ', TS)

    # TODO:
    # Bottle should provide the metadata already,
    # perhaps via the following:
    # from bottle import response, request
    metadata = ast.literal_eval(metadata)

    lipid = metadata['lipid']
    protein = metadata['protein']

    if lipid == "" and protein == "":
        # Starting setup:
        try:
            lipid = BACKEND_DATA['lipids'][0]
            protein = BACKEND_DATA['proteins'][0]
        except:
            print ('Detached EXECUTION')
            print ('This is currently meant for testing only.')
            BACKEND_DATA = independent_execution()
            lipid = BACKEND_DATA['lipids'][0]
            protein = BACKEND_DATA['proteins'][0]


    # For development, let's try to get both the frames and distance_array
    # for a particular lipid selection as an example: 2873
    lipid_id = 2230 # 2873
    gantt_data, categories = get_gantt_app_data(BACKEND_DATA['lipid_contact_frames'], lipid_id)


    # WORKING ON: Table
    table_data = []
    for ix, (lipid_id, freq) in enumerate(BACKEND_DATA['top_lipids']['CHOL']):
        table_data.append({
            "id": ix,
            "lipidID": lipid_id,
            "contactFrequency": freq
        })

    # TODO:
    # Possibly, avoid single point of failure on these dictionary lookups?
    response = {
        "data": {lipid: BACKEND_DATA["data"][protein][lipid]},
        "proteins": BACKEND_DATA['proteins'],
        "lipids": BACKEND_DATA['lipids'],
        "pieData": BACKEND_DATA['pie_data'],
        # "ganttData": BACKEND_DATA['gantt_data'],
        "ganttData": gantt_data,
        # "topLipids": BACKEND_DATA['top_10_lipids'],
        "topLipids": categories,
        "globalTopLipids": BACKEND_DATA['top_lipids'],
        "lipidContactFrames": BACKEND_DATA['lipid_contact_frames'],
        "tableData": table_data
    }
    return response

def start_server(payload=None, debug_bool=False, reloader=True, port=8351):

    # UFCC calls:
    args = payload
    ts = UFCC(args.structure, args.trajectory, add_lipid_types = args.other_lipids)
    ts.contacts.runner.backend = 'prolint_serial'
    ts.contacts.compute(cutoff=args.cutoff)
    payload = ts.contacts.server_payload()

    t, g = sort_lipids(ts)
    payload['top_lipids'] = t
    payload['lipid_contact_frames'] = g

    # Make data accessible globally
    global BACKEND_DATA
    global TS
    BACKEND_DATA = payload
    TS = ts

    debug(debug_bool)
    run(reloader=reloader, host='localhost', port=port)


def independent_execution():
    """
    If we are not calling the server through the ufcc executable, but
    independently, locally for testing purposes, we will load local data file
    and serve that to the dashboard.
    """
    with open(os.path.join(SERVER_PATH, 'girk.json'), 'r') as fp:
        data = json.load(fp)

    pie_data = [{
        "category": "LocalGirk",
        "value": 500,
        "subData": [
            { "category": "CHOL", "value": 300 },
            { "category": "POPE", "value": 150 },
            { "category": "POPS", "value": 50 }
        ]
    }]

    # ganttApp data input requirement
    gantt_data = [{
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
        }
    ]
    top_10_lipids = ['Lipid 1', 'Lipid 2']

    payload = {
        "data": data,
        "proteins": ['LocalGirk'],
        "lipids": list(data['LocalGirk'].keys()),
        "pie_data": pie_data,
        "gantt_data": gantt_data,
        "top_10_lipids": top_10_lipids
    }

    return payload

if __name__ == '__main__':
    start_server(debug_bool=True)