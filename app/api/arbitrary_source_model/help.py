from flask import render_template, url_for, redirect, request, current_app
from .. import api


@api.route('/arbitrary_source_model_help', methods=['GET', 'POST'])
def arbitrary_source_model_help():
    return render_template('arbitrary_source_model/help.html')