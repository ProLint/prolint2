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
    # print (metadata.split(':'))
    if metadata == '1':
        data = dataf('1')
    elif metadata == '2':
        data = dataf('2')
    else:
        data = dataf('3')
    return json.dumps(data)






debug(True)
run(reloader=True, host='localhost', port=8080)

