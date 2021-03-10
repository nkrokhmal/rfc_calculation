import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import io


def show_model(data, dx):
    mat_contents = sio.loadmat(data)
    keys = [key for key in mat_contents.keys() if key[0] != '_']
    pressure_field = mat_contents[keys[0]]
    figure = plt.figure(figsize=(5, 5))
    n = pressure_field.shape[0]
    extent = [- n * dx / 2, n * dx / 2, - n * dx / 2, n * dx / 2]
    im = plt.matshow(np.abs(pressure_field), fignum=figure.number, extent=extent)
    plt.colorbar()
    buffer = io.BytesIO()
    figure.savefig(buffer, format='png')
    buffer.seek(0)
    return buffer