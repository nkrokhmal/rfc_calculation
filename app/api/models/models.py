from flask import render_template, url_for, redirect, session, flash
from .. import api
from ... import db
from ...models import Model


@api.route("/delete_model/<model_id>", methods=['DELETE'])
def delete_model(model_id):
    model = db.session.query(Model)\
        .filter(Model.id == model_id)\
        .first()
    db.session.delete(model)
    db.session.commit()
    flash("Model successfully deleted!", "success")
    return redirect(url_for('.get_models'))


@api.route("/get_models", methods=['GET', 'POST'])
def get_models():
    session.clear()
    data = db.session.query(Model).all()
    return render_template('get_models.html', data=data)
