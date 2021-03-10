from flask import render_template
from . import api


@api.route("/")
def index():
    return render_template("index.html")