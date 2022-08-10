import os
import ast
import json

from bottle import route, run, template, debug, static_file, request

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))

BACKEND_DATA = None
data = None
data_loaded = False

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

@route('/data/:metadata')
def listener(metadata):

    global data_loaded
    global data
    global BACKEND_DATA

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
            BACKEND_DATA = independent_execution()
            lipid = BACKEND_DATA['lipids'][0]
            protein = BACKEND_DATA['proteins'][0]

    response = {
        "data": {lipid: BACKEND_DATA["data"][protein][lipid]},
        "proteins": BACKEND_DATA['proteins'],
        "lipids": BACKEND_DATA['lipids'],
        "pieData": BACKEND_DATA['pie_data']
    }
    return response

def start_server(payload=None, debug_bool=False, reloader=True, port=8351):

    global BACKEND_DATA
    BACKEND_DATA = payload

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
    }

    return payload

if __name__ == '__main__':
    start_server(debug_bool=True)