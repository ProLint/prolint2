import ast
import json
from data import dataf
from bottle import route, run, template, debug, static_file


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
    print ('Receiving request: ')
    print ('metadata', metadata)
    metadata = ast.literal_eval(metadata)

    print (metadata)
    print ("metadata['lipid']", metadata['lipid'])
    print (type(metadata))
    # lipid = metadata.split(':')[0]

    # print ('lipid', lipid)
    lipid = metadata['lipid']
    class PairsHook(dict):
        def __init__(self, pairs):
            key = [x for x in pairs if x[0] == lipid]
            super(PairsHook, self).__init__(key)

    with open('girk.json', 'r') as fp:
        data = json.load(fp, object_pairs_hook=PairsHook)

    # print (metadata.split(':'))
    # if metadata == '1':
    #     data = dataf('1')
    # elif metadata == '2':
    #     data = dataf('2')
    # else:
    #     data = dataf('3')

    # return json.dumps(data)
    print (data.keys())
    return data






debug(True)
run(reloader=True, host='localhost', port=8080)

