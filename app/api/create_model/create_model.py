from flask import render_template, url_for, redirect, request, current_app
from .. import api
from ... import db
from ...models import Model
from ..forms import CreateModelForm
import os
from datetime import datetime
import json
from .focus_beam import build_focused_beam
import scipy.io as sio


@api.route('/create_model', methods=['GET', 'POST'])
def create_model():
    form = CreateModelForm()
    if request.method == 'POST' and form.validate_on_submit():
        model_name = form.model_name.data
        model_path = os.path.join(current_app.config['MODEL_PATH'], f'{model_name}.mat')
        distribution_path = os.path.join(current_app.config['DISTRIBUTION_PATH'], f'{model_name}.png')

        figure, p_flat, z_surf = build_focused_beam(
            r_in=form.radius_of_hole.data,
            r_out=form.radius_of_transducer.data,
            dx=form.spatial_step.data,
            R0=form.curvative_radius.data,
            f=form.frequency.data,
            rho=form.density_of_water.data,
            c=form.speed_of_sound_in_water.data,
            Ampl=form.pressure_amplitude.data)

        figure.savefig(distribution_path)
        sio.savemat(model_path, {'P_ax': p_flat})

        params = {
            "dx": form.spatial_step.data,
            "frequency": form.frequency.data,
            "speed_of_sound": form.speed_of_sound_in_water.data,
            "density_of_medium": form.density_of_water.data,
            "z_surf": z_surf
        }
        model = Model(
            name=model_name,
            path=os.path.join(current_app.config['MODEL_FOLDER'], f'{model_name}.mat'),
            pressure_distribution_path=os.path.join(current_app.config['DISTRIBUTION_FOLDER'], f'{model_name}.png'),
            creation_time=datetime.utcnow(),
            params=json.dumps(params),
            status_id=1
        )
        db.session.add(model)
        db.session.commit()
        return render_template('create_model.html', form=form, figure=model.pressure_distribution_path)
    return render_template('create_model.html', form=form, figure=None)


