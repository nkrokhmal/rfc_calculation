from flask import render_template, url_for, redirect, request, current_app
from .. import api


@api.route('/focused_source_model_help', methods=['GET', 'POST'])
def focused_source_model_help():
    return render_template('focused_source_model/help.html')