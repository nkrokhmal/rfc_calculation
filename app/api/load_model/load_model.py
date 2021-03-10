from flask import request, current_app, render_template
import os
from ...models import Model
from ... import db
from ..forms import UploadModelForm
from .. import api
from .show_model import show_model
from datetime import datetime
import io
import json


@api.route('/load_model', methods=["GET", "POST"])
def load_model():
    form = UploadModelForm()
    if request.method == 'POST' and form.validate_on_submit():
        model_name = form.model_name.data
        model_path = os.path.join(current_app.config['MODEL_PATH'], f'{model_name}.mat')
        distribution_path = os.path.join(current_app.config['DISTRIBUTION_PATH'], f'{model_name}.png')

        params = json.dumps(
            {
                'dx': form.dxvalue.data,
                'frequency': form.frequency.data,
                'speed_of_sound': form.speed_of_sound.data,
                'density_of_medium': form.density_of_medium.data,
                'z_surf': form.z_surf.data,
            }
        )
        file_bytes = io.BytesIO(request.files["input_file"].read())
        figure = show_model(file_bytes, form.dxvalue.data)
        with open(distribution_path, 'wb') as file:
            file.write(figure.getbuffer())

        file_bytes.seek(0)
        with open(model_path, 'wb') as file:
            file.write(file_bytes.getbuffer())

        model = Model(
            name=model_name,
            path=os.path.join(current_app.config['MODEL_FOLDER'], f'{model_name}.mat'),
            pressure_distribution_path=os.path.join(current_app.config['DISTRIBUTION_FOLDER'], f'{model_name}.png'),
            creation_time=datetime.utcnow(),
            params=params,
            status_id=1
        )
        db.session.add(model)
        db.session.commit()
        return render_template('load_model.html', form=form, figure=model.pressure_distribution_path)
    return render_template('load_model.html', form=form, figure=None)
