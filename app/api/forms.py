from flask_wtf import FlaskForm
from flask_wtf.file import FileRequired, FileAllowed
from wtforms import StringField, SubmitField, SelectField, FloatField, FileField
from wtforms.validators import Required, InputRequired

from app import db
from ..models import Model


class NameForm(FlaskForm):
    name = StringField('What is your name?', validators=[Required()])
    curname = StringField('What is your curname?', validators=[Required()])
    submit = SubmitField('Submit')


class ImageForm(FlaskForm):
    url = StringField('What is your avatar?', validators=[Required()])
    submit = SubmitField('Submit')


class UploadModelForm(FlaskForm):
    validators = [
        FileRequired(message='There was no file!'),
        FileAllowed(['mat'], message='Must be a mat file!')
    ]

    input_file = FileField('', validators=validators)
    model_name = StringField("Enter model name:", validators=[Required()])
    dxvalue = FloatField("Enter value dx, m: ", default=0.001)
    frequency = FloatField("Enter value of frequency", default=1000000)
    speed_of_sound = FloatField("Enter value of speed of sound", default=1500.0)
    density_of_medium = FloatField("Enter of density of medium", default=1000.0)
    z_surf = FloatField("Enter Z coordinate of hologram, m", default=0.0)
    submit = SubmitField(label="Submit")


class ScattererForm(FlaskForm):
    model_name = SelectField('Field name', coerce=int, validators=[InputRequired()])
    radius = FloatField("Enter value radius, m", default=0.0001)
    longitudinal = FloatField("Enter longitudinal speed of sound, m/s", default=2620.0)
    transverse = FloatField("Enter transverse speed of sound, m/s", default=1080.0)
    density_of_scatter = FloatField("Enter density of scatterer, kg/m^3", default=1125.0)

    type_value = StringField("Enter type of coordinates (X, Y or Z)", default='Z')
    from_value = FloatField("Enter begin coordinate value", default=-0.02)
    to_value = FloatField("Enter end coordinate value", default=0.02)
    step = FloatField("Enter step value", default=0.001)

    type_value2 = StringField("Enter type of coordinates (X, Y or Z)", default='X')
    from_value2 = FloatField("Enter begin coordinate value", default=-0.02)
    to_value2 = FloatField("Enter end coordinate value", default=0.02)
    step2 = FloatField("Enter step value", default=0.001)

    submit = SubmitField(label="Submit")

    def __init__(self, *args, **kwargs):
        super(ScattererForm, self).__init__(*args, **kwargs)
        self.model_names = db.session.query(Model).all()
        self.model_name.choices = list(enumerate(set([x.name for x in self.model_names])))


class CreateModelForm(FlaskForm):
    model_name = StringField('Enter model name', validators=[Required()])
    radius_of_hole = FloatField('Enter value of radius of hole', default=0.)
    radius_of_transducer = FloatField('Enter value of radius of transducer', default=0.05)
    spatial_step = FloatField('Enter value of dx', default=0.001)
    curvative_radius = FloatField('Enter value of curvative radius', default=0.07)
    frequency = FloatField('Enter value of frequency', default=1000000.)
    density_of_water = FloatField('Enter density of medium', default=1000.)
    speed_of_sound_in_water = FloatField('Enter value of speed of sound', default=1500.)
    pressure_amplitude = FloatField('Enter value of pressure amplitude', default=1.)
    submit = SubmitField(label="Submit")


class ModelResultsForm(FlaskForm):
    model_name = SelectField('Field name', coerce=int, validators=[InputRequired()], default=None)
    submit = SubmitField(label="Search")

    def __init__(self, *args, **kwargs):
        super(ModelResultsForm, self).__init__(*args, **kwargs)
        self.model_names = db.session.query(Model).all()
        self.model_name.choices = list(enumerate(set([x.name for x in self.model_names])))