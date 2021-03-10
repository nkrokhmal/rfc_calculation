from flask import request, current_app, render_template
import os
from ...models import Model, ModelResult
from ... import db
from ..forms import ScattererForm
from .. import api
import json
from ...utils.rfc_client import Object, Wave, Coordinates, Spectrum, Points
import numpy as np
import time


def get_choice_data(f):
    return dict(f.choices).get(f.data)


@api.route("/scatterer", methods=['GET', 'POST'])
def scatterer():
    form = ScattererForm()
    if request.method == 'POST' and form.validate_on_submit():
        cur_time = time.strftime('%Y%m%d%H%M%S')
        model = db.session.query(Model).filter(Model.name == get_choice_data(form.model_name)).first()
        obj = Object(
            a=form.radius.data,
            rho=form.density_of_scatter.data,
            c_l=form.longitudinal.data,
            c_t=form.transverse.data
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
        coordinates_range = np.arange(form.from_value.data, form.to_value.data, form.step.data)
        coordinates_range2 = np.arange(form.from_value2.data, form.to_value2.data, form.step2.data)

        def set_coordinates(key):
            if form.type_value.data == key:
                return coordinates_range
            elif form.type_value2.data == key:
                return coordinates_range2
            else:
                return np.array([0.0])
        coordinates = Coordinates(
            x=set_coordinates('X'),
            y=set_coordinates('Y'),
            z=set_coordinates('Z'),
            z_surf=model_params['z_surf'])

        points = Points(coordinates, obj, wave, spectrum, os.path.join(current_app.config['DATA_FOLDER'], model.path))
        force = points.calculate_force()
        fig1, fig2, fig3 = points.build_rad_force(force)

        fig1.savefig(os.path.join(current_app.config['FORCE_PATH'], f'{model.name}_force1_{cur_time}.png'))
        fig2_path = None
        fig3_path = None
        if fig2:
            fig2_path = os.path.join(current_app.config['FORCE_FOLDER'], f'{model.name}_force2_{cur_time}.png')
            fig2.savefig(os.path.join(current_app.config['FORCE_PATH'], f'{model.name}_force2_{cur_time}.png'))
        if fig3:
            fig3_path = os.path.join(current_app.config['FORCE_FOLDER'], f'{model.name}_force3_{cur_time}.png')
            fig3.savefig(os.path.join(current_app.config['FORCE_PATH'], f'{model.name}_force3_{cur_time}.png'))
        print(fig3_path, fig2_path)
        np.savetxt(os.path.join(current_app.config['FORCE_PATH'], f'{model.name}_force_{cur_time}.txt'), force)

        scatterer_params = {
            'Radius': form.radius.data,
            'LongitudinalSpeed': form.longitudinal.data,
            'TransverseSpeed': form.transverse.data,
            'DensityOfScatterer': form.density_of_scatter.data,
            'Frequency': model_params['frequency'],
            'SpeedOfSound': model_params['speed_of_sound'],
            'DensityOfMedium': model_params['density_of_medium'],
            'Dx': model_params['dx'],
            'Type': form.type_value.data,
            'From': form.from_value.data,
            'To': form.to_value.data,
            'Step': form.step.data,
            'ZSurf': model_params['z_surf']
        }

        model_result = ModelResult(
            # x_force=x_force,
            # y_force=y_force,
            # z_force=z_force,
            force_data_path=os.path.join(current_app.config['FORCE_FOLDER'], f'{model.name}_force_{cur_time}.txt'),
            force_image_path=os.path.join(current_app.config['FORCE_FOLDER'], f'{model.name}_force1_{cur_time}.png'),
            model_id=model.id,
            model_params=json.dumps(scatterer_params),
            status_id=1)
        db.session.add(model_result)
        db.session.commit()
        return render_template('scatterer.html', form=form,
                               figure1=model_result.force_image_path,
                               figure2=fig2_path,
                               figure3=fig3_path,
                               data=model_result.force_data_path)
    return render_template('scatterer.html', form=form, figure=None, data=None)
