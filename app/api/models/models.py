from flask import render_template, url_for, redirect
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
    return redirect(url_for('.get_models'))


@api.route("/get_models", methods=['GET'])
def get_models():
    data = db.session.query(Model).all()
    return render_template('get_models.html', data=data)
