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
    # protein = metadata['protein']

    if not data_loaded:

        class PairsHook(dict):
            def __init__(self, pairs):
                key = [x for x in pairs if x[0] == lipid]
                super(PairsHook, self).__init__(key)

        with open('girk.json', 'r') as fp:
            # data = json.load(fp, object_pairs_hook=PairsHook)
            data = json.load(fp)

        data_loaded = True

    sliced_data = data[lipid]
    return {lipid: sliced_data}



debug(True)
run(reloader=True, host='localhost', port=8080)

