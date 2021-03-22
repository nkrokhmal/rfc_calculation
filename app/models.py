from sqlalchemy import Column, DateTime, String, Integer, func, ForeignKey
from sqlalchemy.orm import backref
from . import db


class Status(db.Model):
    __tablename__ = 'statuses'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String)
    models = db.relationship("Model", backref=backref("status", lazy="subquery"))

    @staticmethod
    def insert_default_status():
        for name in ['Active', 'Disabled']:
            status = Status(name=name)
            db.session.add(status)

        try:
            db.session.commit()
        except Exception as e:
            print(f'Exception occurred {e}')
            db.session.rollback()


class Model(db.Model):
    __tablename__ = 'models'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, unique=True)
    path = db.Column(db.String)
    pressure_distribution_path = db.Column(db.String)
    creation_time = db.Column(db.DateTime, default=func.now())
    params = db.Column(db.String)
    status_id = db.Column(db.Integer, db.ForeignKey('statuses.id'), nullable=True)

    model_results = db.relationship("ModelResult", backref=backref("model", lazy="subquery"))

    def to_json(self):
        return {
            'id': self.id,
            'name': self.model_name,
            'path': self.model_path,
            'pressure_distribution_path': self.pressure_distribution_path,
            'creation_time': self.creation_time,
            'params': self.params,
            'status_id': self.status_id
        }


class ModelResult(db.Model):
    __tablename__ = 'model_results'
    id = db.Column(db.Integer, primary_key=True)
    x_force = db.Column(db.Float)
    y_force = db.Column(db.Float)
    z_force = db.Column(db.Float)
    force_data_path = db.Column(db.String)
    force_image = db.Column(db.String)
    model_params = db.Column(db.String)
    status_id = db.Column(db.Integer, db.ForeignKey('statuses.id'))
    model_id = db.Column(db.Integer, db.ForeignKey('models.id'))

    def to_json(self):
        return {
            'id': self.id,
            'x_force': self.x_force,
            'y_force': self.y_force,
            'z_force': self.x_force,
            'force_data_path': self.force_data_path,
            'force_image_path': self.force_image_path,
            'model_params': self.model_params,
            'status_id': self.status_id,
        }
