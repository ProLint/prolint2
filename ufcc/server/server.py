import ast
import json
from bottle import route, run, template, debug, static_file

# rendered_data = {
#     "proteins": {
#         "name": None,
#         "lipids": []
#     },
# }
data = None
data_loaded = False

@route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root='static')

@route('/')
def index():
    return template('home.tpl')

@route('/app')
def app():
    return static_file('index.html', root='.')

@route('/data/:metadata')
def listener(metadata):

    global data_loaded
    global data

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

        with open('girk.json', 'r') as fp:
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



debug(True)
run(reloader=True, host='localhost', port=8351)

