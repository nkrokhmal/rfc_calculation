from flask import request, url_for, render_template, redirect, flash, session
from ...models import Scatterer
from ... import db
from ..forms import SphericalScattererParametersForm
from .. import api


@api.route("/get_scatterer_params", methods=['GET', 'POST'])
def get_scatterer_params():
    session.clear()
    scatterers = db.session.query(Scatterer).all()
    return render_template('scatterer_params/get_scatterer_params.html', scatterers=scatterers)


@api.route("/add_scatterer_params", methods=['GET', 'POST'])
def add_scatterer_params():
    form = SphericalScattererParametersForm()
    if form.validate_on_submit():
        scatter = Scatterer(
            name=form.name.data,
            longitudinal_velocity=form.longitudinal_velocity.data,
            shear_velocity=form.shear_velocity.data,
            density=form.density.data,
        )
        db.session.add(scatter)
        db.session.commit()

        return redirect(url_for(".get_scatterer_params"))

    return render_template('scatterer_params/add_scatterer_params.html', form=form)


@api.route("/edit_scatterer_params/<int:scatterer_id>", methods=['GET', 'POST'])
def edit_scatterer_params(scatterer_id):
    form = SphericalScattererParametersForm()
    scatterer = db.session.query(Scatterer).get_or_404(scatterer_id)
    if request.method == 'POST' and form.validate_on_submit():
        scatterer.name = form.name.data
        scatterer.longitudinal_velocity = form.longitudinal_velocity.data
        scatterer.shear_velocity = form.shear_velocity.data
        scatterer.density = form.density.data
        db.session.commit()
        return redirect(url_for(".get_scatterer_params"))

    form.name.data = scatterer.name
    form.longitudinal_velocity.data = scatterer.longitudinal_velocity
    form.shear_velocity.data = scatterer.shear_velocity
    form.density.data = scatterer.density

    return render_template('scatterer_params/edit_scatterer_params.html', form=form, scatterer_id=scatterer.id)


@api.route("/delete_scatterer_params/<int:scatterer_id>", methods=["DELETE"])
def delete_scatterer_params(scatterer_id):
    scatterer = db.session.query(Scatterer).get_or_404(scatterer_id)
    if scatterer:
        db.session.delete(scatterer)
        db.session.commit()
        flash("Scatterer successfully deleted!", "success")
    return redirect(url_for(".get_scatterer_params"))