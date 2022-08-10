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

    metadata = ast.literal_eval(metadata)

    lipid = metadata['lipid']
    protein = metadata['protein']

    if lipid == "" and protein == "":
        # Starting setup:
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

if __name__ == '__main__':
    start_server(debug_bool=True)