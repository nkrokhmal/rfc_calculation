import os
from flask_migrate import MigrateCommand, Migrate
from flaskwebgui import FlaskUI

if os.path.exists('.env'):
    for line in open('.env'):
        var = line.strip().split('=')
        if len(var) == 2:
            os.environ[var[0]] = var[1]

from flask_script import Manager
from app import create_app

app, db = create_app(os.getenv('FLASK_CONFIG') or 'default')
manager = Manager(app)
migrate = Migrate(app, db, render_as_batch=True)
manager.add_command("db", MigrateCommand)
# ui = FlaskUI(app)
#         admin.add_view(ModelView(obj, db.session))

if __name__ == '__main__':
    manager.run()
    # app.run(debug=True, port=5000, threaded=True, host='0.0.0.0')