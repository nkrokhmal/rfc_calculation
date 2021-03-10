from . import api
from flask import send_from_directory


@api.route('/data/<path:path>')
def send_data(path):
    return send_from_directory('data', path)