from flask_wtf import FlaskForm
from flask_wtf.file import FileRequired, FileAllowed
from wtforms import StringField, SubmitField, SelectField, FloatField, FileField, ValidationError
from wtforms.validators import Required, InputRequired

from app import db
from ..models import Model, Scatterer


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
    material = SelectField('Material name', coerce=int, validators=[InputRequired()])
    radius = FloatField("Enter value radius, m", default=0.0001)

    from_value_x = FloatField("Enter begin coordinate value", default=-0.02)
    to_value_x = FloatField("Enter end coordinate value", default=0.02)
    step_x = FloatField("Enter step value", default=0.001)

    from_value_y = FloatField("Enter begin coordinate value", default=-0.02)
    to_value_y = FloatField("Enter end coordinate value", default=0.02)
    step_y = FloatField("Enter step value", default=0.001)

    from_value_z = FloatField("Enter begin coordinate value", default=-0.02)
    to_value_z = FloatField("Enter end coordinate value", default=0.02)
    step_z = FloatField("Enter step value", default=0.001)

    submit = SubmitField(label="Submit")

    def validate_to_value_x(self, field):
        if field.data < self.from_value_x.data:
            raise ValidationError("End coordinate value must be greater than begin coordinate value!")

    def validate_to_value_y(self, field):
        if field.data < self.from_value_y.data:
            raise ValidationError("End coordinate value must be greater than begin coordinate value!")

    def validate_to_value_z(self, field):
        if field.data < self.from_value_z.data:
            raise ValidationError("End coordinate value must be greater than begin coordinate value!")

    def validate_step_x(self, field):
        if self.from_value_x.data < self.to_value_x.data:
            if field.data > self.to_value_x.data - self.from_value_x.data:
                raise ValidationError("Step greater than difference between begin and end coordinates")

    def validate_step_y(self, field):
        if self.from_value_y.data < self.to_value_y.data:
            if field.data > self.to_value_y.data - self.from_value_y.data:
                raise ValidationError("Step greater than difference between begin and end coordinates")

    def validate_step_z(self, field):
        if self.from_value_z.data < self.to_value_z.data:
            if field.data > self.to_value_z.data - self.from_value_z.data:
                raise ValidationError("Step greater than difference between begin and end coordinates")

    def __init__(self, *args, **kwargs):
        super(ScattererForm, self).__init__(*args, **kwargs)
        self.model_names = db.session.query(Model).all()
        self.model_name.choices = list(enumerate(set([x.name for x in self.model_names])))
        self.materials = db.session.query(Scatterer).all()
        self.material.choices = list(enumerate(set([x.name for x in self.materials])))


# todo: rename
class CreateModelForm(FlaskForm):
    model_name = StringField('Enter model name', validators=[Required()])
    radius_of_hole = FloatField('Enter value of radius of hole', default=0.)
    radius_of_transducer = FloatField('Enter value of radius of transducer', default=50)
    spatial_step = FloatField('Enter value of dx', default=1)
    curvative_radius = FloatField('Enter value of curvative radius', default=70)
    frequency = FloatField('Enter value of frequency', default=1.)
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


class SphericalScattererParametersForm(FlaskForm):
    name = StringField('Field name', validators=[Required()])
    longitudinal_velocity = FloatField('Longitudinal velocity, m/s', validators=[Required()])
    shear_velocity = FloatField('Shear velocity, m/s', validators=[Required()])
    density = FloatField('Density, kg/m^3', validators=[Required()])
