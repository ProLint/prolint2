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

    # TODO:
    # Bottle should provide the metadata already,
    # perhaps via the following:
    # from bottle import response, request
    metadata = ast.literal_eval(metadata)

    lipid = metadata['lipid']
    protein = metadata['protein']

    if lipid == "" and protein == "":
        # Starting setup:
        lipid = "CHOL"
        protein = "GIRK"

    if not data_loaded:

        # TODO:
        # determine if the PairsHook class is/will be needed
        # class PairsHook(dict):
        #     def __init__(self, pairs):
        #         key = [x for x in pairs if x[0] == lipid]
        #         super(PairsHook, self).__init__(key)

        with open('girk.json', 'r') as fp:
            # data = json.load(fp, object_pairs_hook=PairsHook)
            data = json.load(fp)

        data_loaded = True

    # TODO:
    # get the protein from the server
    sliced_data = data['Protein0'][lipid]

    # PieChart App Input Data:
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
    }, # => the section below shows how subsequent proteins should be included.
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
    # This is the object we will send to the front end.
    # It should have everything the apps need to work.
    # "data" is required by the radarApp
    # "proteins", "lipids", and "pieData" are needed by the pieApp # TODO: maybe we should simplify these requirements?
    #
    response = {
        "data": {lipid: sliced_data},
        "proteins": ['Protein0'],
        "lipids": list(data['Protein0'].keys()),
        "pieData": pie_data,
        "ganttData": gantt_data
    }
    return response
    # return {lipid: sliced_data}



debug(True)
run(reloader=True, host='localhost', port=8351)

