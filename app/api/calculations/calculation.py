from flask import request, current_app, render_template
import os
from ...models import Model, ModelResult, Scatterer
from ... import db
from ..forms import ScattererForm
from .. import api
import json
from ...utils.rfc_client import Object, Wave, Coordinates, Spectrum, Points
import numpy as np
import time
import plotly


def get_choice_data(f):
    return dict(f.choices).get(f.data)


def create_arrange_with_endpoint(beg, end, step):
    return np.arange(beg, end + step, step)


@api.route("/calculation", methods=['GET', 'POST'])
def scatterer():
    form = ScattererForm()
    if request.method == 'POST' and form.validate_on_submit():
        cur_time = time.strftime('%Y%m%d%H%M%S')
        model = db.session.query(Model).filter(Model.name == get_choice_data(form.model_name)).first()
        material = db.session.query(Scatterer).filter(Scatterer.name == get_choice_data(form.material)).first()
        obj = Object(
            a=form.radius.data,
            rho=material.density,
            c_l=material.longitudinal_velocity,
            c_t=material.shear_velocity
        )
        model_params = json.loads(model.params)
        wave = Wave(
            f=model_params['frequency'],
            c=model_params['speed_of_sound'],
            rho=model_params['density_of_medium'],
        )
        spectrum = Spectrum(
            dx=model_params['dx']
        )
        coordinates_range_x = create_arrange_with_endpoint(form.from_value_x.data, form.to_value_x.data, form.step_x.data) if \
            form.to_value_x.data > form.from_value_x.data else np.array([form.from_value_x.data])
        coordinates_range_y = create_arrange_with_endpoint(form.from_value_y.data, form.to_value_y.data, form.step_y.data) if \
            form.to_value_y.data > form.from_value_y.data else np.array([form.from_value_y.data])
        coordinates_range_z = create_arrange_with_endpoint(form.from_value_z.data, form.to_value_z.data, form.step_z.data) if \
            form.to_value_z.data > form.from_value_z.data else np.array([form.from_value_z.data])

        coordinates = Coordinates(
            x=coordinates_range_x,
            y=coordinates_range_y,
            z=coordinates_range_z,
            z_surf=model_params['z_surf'])

        points = Points(coordinates, obj, wave, spectrum, os.path.join(current_app.config['DATA_FOLDER'], model.path))
        force = points.calculate_force()
        figure = points.build_rad_force(force)
        figure_params = json.dumps(figure, cls=plotly.utils.PlotlyJSONEncoder)
        np.savetxt(os.path.join(current_app.config['FORCE_PATH'], f'{model.name}_force_{cur_time}.txt'), force)

        scatterer_params = {
            'Radius': form.radius.data,
            'LongitudinalSpeed': material.longitudinal_velocity,
            'TransverseSpeed': material.shear_velocity,
            'DensityOfScatterer': material.density,
            'Frequency': model_params['frequency'],
            'SpeedOfSound': model_params['speed_of_sound'],
            'DensityOfMedium': model_params['density_of_medium'],
            'Dx': model_params['dx'],
            'From_x': form.from_value_x.data,
            'To_x': form.to_value_x.data,
            'Step_x': form.step_x.data,
            'From_y': form.from_value_y.data,
            'To_y': form.to_value_y.data,
            'Step_y': form.step_y.data,
            'From_z': form.from_value_z.data,
            'To_z': form.to_value_z.data,
            'Step_z': form.step_z.data,
            'ZSurf': model_params['z_surf']
        }

        model_result = ModelResult(
            force_data_path=os.path.join(current_app.config['FORCE_FOLDER'], f'{model.name}_force_{cur_time}.txt'),
            force_image=figure_params,
            model_id=model.id,
            model_params=json.dumps(scatterer_params),
            status_id=1)
        db.session.add(model_result)
        db.session.commit()
        return render_template('calculations/calculation.html', form=form,
                               figure=figure_params)
    return render_template('calculations/calculation.html', form=form, figure=None)
