from flask import render_template, url_for, redirect, request
from .. import api
from ... import db
from ...models import ModelResult, Model
from ..forms import ModelResultsForm


@api.route("/delete_model_result/<int:model_result_id>", methods=['DELETE'])
def delete_model_result(model_result_id):
    model_result = db.session.query(ModelResult) \
        .filter(ModelResult.id == model_result_id) \
        .first()
    db.session.delete(model_result)
    db.session.commit()
    return redirect(url_for('.get_models'))


@api.route("/get_model_result/<int:id>", methods=['GET', 'POST'])
def get_model_result(id):
    model = db.session.query(Model).filter(Model.id == id).first()
    model_results = db.session.query(ModelResult)\
        .filter(ModelResult.model_id == id)\
        .all()
    return render_template('get_model_result.html', model_results=model_results, model=model)
