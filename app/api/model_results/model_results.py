from flask import render_template, url_for, redirect, request
from .. import api
from ... import db
from ...models import ModelResult
from ..forms import ModelResultsForm


@api.route("/delete_model_result/<model_result_id>", methods=['DELETE'])
def delete_model_result(model_result_id):
    model_result = db.session.query(ModelResult) \
        .filter(ModelResult.id == model_result_id) \
        .first()
    db.session.remove(model_result)
    db.session.commit()
    return redirect(url_for('.get_model_result'))


@api.route("/get_model_result", methods=['GET', 'POST'])
def get_model_result():
    form = ModelResultsForm()
    if form.validate_on_submit():
        model_results = db.session.query(ModelResult)\
            .filter(ModelResult.model.name == form.model_name.data)\
            .all()
        return render_template('get_model_result.html', form=form, model_results=model_results)
    model_results = db.session.query(ModelResult).all()
    return render_template('get_model_result.html', form=form, model_results=model_results)