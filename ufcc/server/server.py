import os
import ast
import json
from bottle import route, run, template, debug, static_file, request

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))
BACKEND_DATA = None
# rendered_data = {
#     "proteins": {
#         "name": None,
#         "lipids": []
#     },
# }
# data = None
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
    global BACKEND_DATA
    print ('BACKEND_DATA', BACKEND_DATA)

    metadata = ast.literal_eval(metadata)

    lipid = metadata['lipid']
    protein = metadata['protein']

    if lipid == "" and protein == "":
        # Starting setup:
        lipid = "CHOL"
        protein = "GIRK"

    if not data_loaded:

        # class PairsHook(dict):
        #     def __init__(self, pairs):
        #         key = [x for x in pairs if x[0] == lipid]
        #         super(PairsHook, self).__init__(key)

        with open(os.path.join(SERVER_PATH, 'girk.json'), 'r') as fp:
            # data = json.load(fp, object_pairs_hook=PairsHook)
            data = json.load(fp)

        data_loaded = True

    sliced_data = data['Protein0'][lipid]

    # Value can be system data: e.g. the ratio of the different lipids, but in that case all
    # values for all different proteins would be the same (not necessarily a bad thing)
    # Values can also be relative ratio of contacts with the different lipids, in which case
    # they are specific for any protein (then again, without normalization it's unclear how
    # useful this information is.)
    pie_data = [{
        "category": "Protein0",
        "value": 500,
        "subData": [
            { "category": "CHOL", "value": 300 },
            { "category": "POPE", "value": 150 },
            { "category": "POPS", "value": 50 }
        ]
    },
    # {
    #     "category": "Protein1",
    #     "value": 300,
    #     "subData": [
    #         { "category": "CHOL", "value": 100 },
    #         { "category": "POPE", "value": 150 },
    #         { "category": "POPS", "value": 50 }

    #     ]
    # }
    ]

    response = {
        "data": {lipid: sliced_data},
        "proteins": ['Protein0'],
        "lipids": list(data['Protein0'].keys()),
        "pieData": pie_data
    }
    return response
    # return {lipid: sliced_data}

def start_server(payload=None, debug_bool=False, reloader=True, port=8351):

    # import requests
    global BACKEND_DATA
    BACKEND_DATA = payload
    # global data
    # data = payload
    # response = requests.request("POST", url, data=payload, headers=headers)

    debug(debug_bool)
    run(reloader=reloader, host='localhost', port=port)

if __name__ == '__main__':
    start_server(debug_bool=True)