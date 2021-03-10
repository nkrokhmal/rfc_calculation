from flask import Blueprint

api = Blueprint('api', __name__)

from . import urls
from . import index
from .create_model import *
from .model_results import *
from .models import *
from .load_model import *
from .scatterer import *
