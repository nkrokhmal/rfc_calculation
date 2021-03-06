"""empty message

Revision ID: a160047fbeeb
Revises: 
Create Date: 2021-06-03 12:12:48.280267

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'a160047fbeeb'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('scatterers',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(), nullable=True),
    sa.Column('longitudinal_velocity', sa.Float(), nullable=True),
    sa.Column('shear_velocity', sa.Float(), nullable=True),
    sa.Column('density', sa.Float(), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('statuses',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('models',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(), nullable=True),
    sa.Column('path', sa.String(), nullable=True),
    sa.Column('pressure_distribution_path', sa.String(), nullable=True),
    sa.Column('creation_time', sa.DateTime(), nullable=True),
    sa.Column('params', sa.String(), nullable=True),
    sa.Column('status_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['status_id'], ['statuses.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('name')
    )
    op.create_table('model_results',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('x_force', sa.Float(), nullable=True),
    sa.Column('y_force', sa.Float(), nullable=True),
    sa.Column('z_force', sa.Float(), nullable=True),
    sa.Column('force_data_path', sa.String(), nullable=True),
    sa.Column('force_image', sa.String(), nullable=True),
    sa.Column('model_params', sa.String(), nullable=True),
    sa.Column('status_id', sa.Integer(), nullable=True),
    sa.Column('model_id', sa.Integer(), nullable=True),
    sa.ForeignKeyConstraint(['model_id'], ['models.id'], ),
    sa.ForeignKeyConstraint(['status_id'], ['statuses.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('model_results')
    op.drop_table('models')
    op.drop_table('statuses')
    op.drop_table('scatterers')
    # ### end Alembic commands ###
