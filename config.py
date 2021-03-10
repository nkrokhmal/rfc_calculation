import os

basedir = os.path.abspath(os.path.dirname(__file__))
SQLITE_PATH = os.path.join(basedir, "data.sqlite")


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY') or "Alisa"
    MODEL_PATH = 'app/data/models'
    MODEL_FOLDER = 'models'
    DISTRIBUTION_PATH = 'app/data/distributions'
    DISTRIBUTION_FOLDER = 'distributions'
    TMP_PATH = 'app/data/tmp'
    FORCE_PATH = 'app/data/force'
    FORCE_FOLDER = 'force'

    DATA_FOLDER = 'app/data'

    @staticmethod
    def init_app(app):
        pass


class DebugConfig(Config):
    SQL_CONNECTION_STRING = (os.environ.get('SQLALCHEMY_DATABASE_URI') or "sqlite:///" + SQLITE_PATH) \
                            + "?check_same_thread=False"


class ProductionConfig(Config):
    SSL_DISABLE = False
    SQLALCHEMY_DATABASE_URI = (os.environ.get('SQLALCHEMY_DATABASE_URI') or "sqlite:///" + SQLITE_PATH) \
                            + "?check_same_thread=False"


config = {
    'debug': DebugConfig,
    'production': ProductionConfig,
    'default': ProductionConfig
}
