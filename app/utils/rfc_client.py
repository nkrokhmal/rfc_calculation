import math
import scipy
import scipy.special as special
from os.path import join as pjoin
import scipy.io as sio
from scipy.fftpack import fftn
import plotly.graph_objs as go
import numpy as np
from plotly.subplots import make_subplots
from scipy.special import factorial


class Object:
    def __init__(self, a=0.0005, rho=1125., c_l=2620., c_t=1080., k_l=0., k_t=0.):
        self.a = a
        self.rho = rho
        self.c_l = c_l
        self.c_t = c_t
        self.sigma = (c_l ** 2 / 2 - c_t ** 2) / (c_l ** 2 - c_t ** 2)
        self.k_l = k_l
        self.k_t = k_t


class Wave:
    def __init__(self, f=2.0e6, c=1500., rho=1000.):
        self.f = f
        self.c = c
        self.k = 2 * math.pi * f / c
        self.rho = rho


class Spectrum:
    def __init__(self, dx=0.00015, dk=0, Nx=0, Ny=0):
        self.dx = dx
        self.dk = dk
        self.Nx = Nx
        self.Ny = Ny


class Coordinates:
    def __init__(self, z=np.arange(-0.01, 0.0105, 0.0005), y=np.arange(-0.005, 0.0055, 0.0001), x=np.array([0.0]), z_surf=None):
        self.x = x
        self.y = y
        self.z = z
        self.N_points = x.shape[0] * y.shape[0] * z.shape[0]
        self.nx = x.shape[0]
        self.ny = y.shape[0]
        self.nz = z.shape[0]
        self.z_serf = z_surf


class Points:
    def __init__(self, coordinates, obj, wave, spectrum, path):
        self.coordinates = coordinates
        self.wave = wave
        self.obj = obj
        self.spectrum = spectrum
        self.obj.k_l = 2 * math.pi * self.wave.f / self.obj.c_l
        self.obj.k_t = 2 * math.pi * self.wave.f / self.obj.c_t
        self.n_global = 5 + int(np.ceil(self.wave.k * self.obj.a));  # 3 * int(np.ceil(self.wave.k * self.obj.a));
        self.n = np.arange(self.n_global + 1)
        self.path = path
        self.points = self.init_points()
        self.init_dif_bessel_kla()
        self.init_dif_bessel_kta()
        self.init_dif_bessel_ka()
        self.init_bessel()
        self.init_c_n()
        self.init_angle_spectrum()
        self.init_wave_number_space()
        self.init_spherical_legandre()

    def init_points(self):
        points = np.array(
            [[x, y, z] for x in self.coordinates.x for y in self.coordinates.y for z in self.coordinates.z])
        return points

    def init_dif_bessel_kla(self):
        self.dif_bessel_kla = self.n * (0. + 1j * 0.)
        self.dif2_bessel_kla = self.n * (0. + 1j * 0.)

    def init_dif_bessel_kta(self):
        self.dif_bessel_kta = self.n * (0. + 1j * 0.)
        self.dif2_bessel_kta = self.n * (0. + 1j * 0.)

    def init_dif_bessel_ka(self):
        self.dif_bessel_ka = self.n * (0. + 1j * 0.)
        self.dif2_bessel_ka = self.n * (0. + 1j * 0.)

    def init_bessel(self):
        self.sph_bessel_kla = scipy.special.spherical_jn(self.n, self.obj.k_l * self.obj.a)
        self.sph_bessel_kta = scipy.special.spherical_jn(self.n, self.obj.k_t * self.obj.a)
        self.sph_bessel_ka = scipy.special.spherical_jn(self.n, self.wave.k * self.obj.a)
        self.sph_hankel = scipy.special.spherical_jn(self.n,
                                                     self.wave.k * self.obj.a) + 1j * scipy.special.spherical_yn(self.n,
                                                                                                                 self.wave.k * self.obj.a)

        self.dif_hankel = self.n * (0. + 1j * 0.)

        self.dif_bessel_kla[0] = - self.sph_bessel_kla[1]
        self.dif_bessel_kta[0] = - self.sph_bessel_kta[1]
        self.dif_bessel_ka[0] = - self.sph_bessel_ka[1]
        self.dif_hankel[0] = - self.sph_hankel[1]

        for k in range(1, self.n_global + 1):
            self.dif_bessel_kla[k] = self.sph_bessel_kla[k - 1] - (k + 1) / (self.obj.k_l * self.obj.a) * \
                                     self.sph_bessel_kla[k]
            if self.obj.k_t == 0:
                self.dif_bessel_kta[k] = self.sph_bessel_kta[k - 1]
            else:
                self.dif_bessel_kta[k] = self.sph_bessel_kta[k - 1] - (k + 1) / (self.obj.k_t * self.obj.a) * \
                                         self.sph_bessel_kta[k]
            self.dif_bessel_ka[k] = self.sph_bessel_ka[k - 1] - (k + 1) / (self.wave.k * self.obj.a) * \
                                    self.sph_bessel_ka[k]
            self.dif_hankel[k] = self.sph_hankel[k - 1] - (k + 1) / (self.wave.k * self.obj.a) * self.sph_hankel[k]

        for i in range(0, self.n_global):
            self.dif2_bessel_kla[i] = i / (2 * i + 1) * self.dif_bessel_kla[i - 1] - (i + 1) / (2 * i + 1) * \
                                      self.dif_bessel_kla[i + 1]
            self.dif2_bessel_kta[i] = i / (2 * i + 1) * self.dif_bessel_kta[i - 1] - (i + 1) / (2 * i + 1) * \
                                      self.dif_bessel_kta[i + 1]
            self.dif2_bessel_ka[i] = i / (2 * i + 1) * self.dif_bessel_ka[i - 1] - (i + 1) / (2 * i + 1) * \
                                     self.dif_bessel_ka[i + 1]

    def init_c_n(self):
        alf = self.sph_bessel_kla - self.obj.k_l * self.obj.a * self.dif_bessel_kla
        bet = (self.n ** 2 + self.n - 2) * self.sph_bessel_kta + (self.obj.k_t * self.obj.a) ** 2 * self.dif2_bessel_kta
        delta = 2 * self.n * (self.n + 1) * self.sph_bessel_kta
        ksi = self.obj.k_l * self.obj.a * self.dif_bessel_kla
        nu = 2 * self.n * (self.n + 1) * (self.sph_bessel_kta - self.obj.k_t * self.obj.a * self.dif_bessel_kta)
        eps = self.obj.k_l ** 2 * self.obj.a ** 2 * (
                    self.sph_bessel_kla * self.obj.sigma / (1 - 2 * self.obj.sigma) - self.dif2_bessel_kla)

        G = self.wave.rho * self.obj.k_t ** 2 * self.obj.a ** 2 / 2 / self.obj.rho * (alf * delta + bet * ksi) / (
                    alf * nu + bet * eps)
        self.c_n = - (G * self.sph_bessel_ka - self.wave.k * self.obj.a * self.dif_bessel_ka) / (
                    G * self.sph_hankel - self.wave.k * self.obj.a * self.dif_hankel)

    def load_file(self):
        mat_fname = pjoin(self.path)
        mat_contents = sio.loadmat(mat_fname)
        return mat_contents

    '''to do: spectrum'''

    def init_angle_spectrum(self):
        mat_contents = self.load_file()

        #         keys = sorted(mat_contents.keys())
        keys = [key for key in mat_contents.keys() if key[0] != '_']
        pressure_field = mat_contents[keys[0]]
        angle_spectrum_0 = scipy.fftpack.fftshift(
            fftn(pressure_field, pressure_field.shape))  ###########################
        self.spectrum.Nx = angle_spectrum_0.shape[0]
        self.spectrum.Ny = angle_spectrum_0.shape[1]
        lin_l = (-1) ** np.arange(self.spectrum.Nx)
        lin_l = lin_l[:, np.newaxis]
        lin_m = lin_l.T
        lin_lm = lin_l @ lin_m
        self.angle_spectrum = angle_spectrum_0.conj() * lin_lm * (4 * np.pi ** 2) / self.spectrum.Nx ** 2

    '''to do: split arrays and angles'''

    def init_wave_number_space(self):
        if (-1) ** self.spectrum.Nx > 0:
            x_array = np.concatenate([np.arange(- self.spectrum.Nx / 2, self.spectrum.Nx / 2, 1)])
        elif (-1) ** self.spectrum.Nx < 0:
            x_array = np.concatenate([np.arange(- (self.spectrum.Nx - 1.) / 2, (self.spectrum.Nx - 1.) / 2 + 1, 1)])
        x_array = x_array[:, np.newaxis]
        x_array = self.spectrum.dx * x_array.astype(float)
        y_array = x_array.copy()
        y_array = y_array.T
        self.r_array = np.sqrt(x_array ** 2 + y_array ** 2)
        #         self.r_array[0,0] = 1e-12

        self.spectrum.dk = 2 * math.pi / (self.spectrum.dx * (self.spectrum.Nx - 1))
        self.kx_array = self.spectrum.dk / self.spectrum.dx * x_array.copy() - 0 * self.spectrum.dk
        self.ky_array = self.spectrum.dk / self.spectrum.dx * y_array.copy() - 0 * self.spectrum.dk
        self.kr_array = np.sqrt(self.kx_array ** 2 + self.ky_array ** 2)
        self.k_window = 0.5 * (np.sign(self.wave.k - self.kr_array) + 1)

        kr2 = np.float_power(self.kr_array, 2)
        kr22 = (self.wave.k ** 2 - kr2) * self.k_window
        self.kz_array = np.float_power(kr22, 1 / 2)

        self.phi_k = (np.arctan2(self.kx_array, self.ky_array)) * self.k_window
        self.cos_th_k = np.float_power((1 - self.kr_array ** 2 / self.wave.k ** 2) * self.k_window, 0.5)
        self.th_k = np.arccos(self.cos_th_k) * self.k_window

        self.phi = (np.arctan2(x_array, y_array))
        self.cos_th = 0 * self.r_array
        self.th = np.arccos(self.cos_th)
        self.power = np.sum(np.sum(1 / 8 / np.pi ** 2 / self.wave.rho / self.wave.c * self.cos_th_k * np.abs(
            self.angle_spectrum) ** 2)) * self.spectrum.dk ** 2

    def calculate_spharm(self, n, m, Y_sph):
        if m >= 0:
            SP = Y_sph[n, m]
        elif m < 0:
            SP = np.conj(Y_sph[n, -m]) * (- 1) ** (- m)
        return SP

    ''' to do: вынести цикл в функцию'''

    def init_spherical_legandre(self):
        Nx = np.size(self.th_k, 0)
        Ny = np.size(self.th_k, 1)

        self.H_nm = np.zeros((self.n_global + 1, 2 * self.n_global + 1, self.coordinates.N_points)) * (0. + 1j * 0.)
        P = np.zeros((Nx, Ny)) * (0. + 1j * 0)
        Y_sph = P.copy()
        Kn = P.copy()

        for nn in range(0, self.n_global + 1):
            if (nn == 0):
                mm = 0
                K = ((2 * nn + 1) / 4 / math.pi * factorial(nn - mm, exact=True) / factorial(nn + mm, exact=True)) ** (
                            1 / 2)
                Kn = K * np.exp(1j * (mm * self.phi_k))
                P = scipy.special.lpmv(mm, nn, self.cos_th_k)
                Y_sph = P * Kn
                for glob_n in range(self.coordinates.N_points):
                    phase_mult = np.exp(1j * self.kx_array * self.points[glob_n, 0] + 1j * self.ky_array * self.points[
                        glob_n, 1] + 1j * self.kz_array * self.points[glob_n, 2])
                    s_newpoint = phase_mult * self.angle_spectrum * self.k_window
                    tmp = Y_sph.conj() * s_newpoint * self.spectrum.dk ** 2
                    self.H_nm[0, 0, glob_n] = tmp.sum(axis=1).sum()
            else:
                for mm in range(0, nn + 1):
                    K = ((2 * nn + 1) / 4 / math.pi * factorial(nn - mm, exact=True) / factorial(nn + mm,
                                                                                                 exact=True)) ** (1 / 2)
                    Kn = K * np.exp(1j * (mm * self.phi_k))
                    P = scipy.special.lpmv(mm, nn, self.cos_th_k)
                    Y_sph_plus = P * Kn
                    Y_sph_minus = (P * Kn).conj() * (-1) ** mm
                    for glob_n in range(self.coordinates.N_points):
                        if mm == 0:
                            phase_mult = np.exp(
                                1j * self.kx_array * self.points[glob_n, 0] + 1j * self.ky_array * self.points[
                                    glob_n, 1] + 1j * self.kz_array * self.points[glob_n, 2])
                            s_newpoint = phase_mult * self.angle_spectrum * self.k_window
                            tmp_0 = Y_sph_plus.conj() * s_newpoint * self.spectrum.dk ** 2
                            self.H_nm[nn, 0, glob_n] = tmp_0.sum(axis=1).sum()
                        else:
                            phase_mult = np.exp(
                                1j * self.kx_array * self.points[glob_n, 0] + 1j * self.ky_array * self.points[
                                    glob_n, 1] + 1j * self.kz_array * self.points[glob_n, 2])
                            s_newpoint = phase_mult * self.angle_spectrum * self.k_window
                            tmp_minus = Y_sph_minus.conj() * s_newpoint * self.spectrum.dk ** 2
                            tmp_plus = Y_sph_plus.conj() * s_newpoint * self.spectrum.dk ** 2
                            self.H_nm[nn, -mm, glob_n] = tmp_minus.sum(axis=1).sum()
                            self.H_nm[nn, mm, glob_n] = tmp_plus.sum(axis=1).sum()

    def calculate_force(self):
        force = np.zeros((self.coordinates.N_points, 3)) * 1.0
        accuracy_z = np.zeros((self.coordinates.N_points, 1)) * 1.0
        Nx = np.size(self.th_k, 0)
        Ny = np.size(self.th_k, 1)
        #         scat_p = np.zeros((self.coordinates.N_points, Nx, Ny)) * (1 + 1j)
        #         p_newpoint = scat_p.copy()
        #         H_nm = np.zeros((self.n_global + 1, 2 * self.n_global + 1, 1)) * (0. + 1j * 0.)

        for glob_n in range(self.coordinates.N_points):  #
            f_nx = np.zeros((self.n_global, 1)) * (0. + 1j * 0.)
            f_ny = np.zeros((self.n_global, 1)) * (0. + 1j * 0.)
            f_nz = np.zeros((self.n_global, 1)) * (0. + 1j * 0.)
            force_x = 0.
            force_y = 0.
            force_z = 0.
            for nn in range(self.n_global):
                psi = (1 + 2 * self.c_n[nn]) * (1 + 2 * self.c_n[nn + 1].conj()) - 1
                f_x = np.zeros((2 * self.n_global + 1, 1)) * (0. + 1j * 0.)
                f_y = f_x.copy()
                f_z = f_x.copy()
                for mm in range(-nn, nn + 1):
                    A_nm = np.sqrt((nn + mm + 1) * (nn + mm + 2) / (2 * nn + 1) / (2 * nn + 3))
                    B_nm = np.sqrt((nn + mm + 1) * (nn - mm + 1) / (2 * nn + 1) / (2 * nn + 3))
                    # print(B_nm)
                    f_x[mm + nn] = A_nm * (
                                self.H_nm[nn, mm, glob_n] * self.H_nm[nn + 1, mm + 1, glob_n].conj() - self.H_nm[
                            nn, - mm, glob_n] * self.H_nm[nn + 1, - mm - 1, glob_n].conj())

                    f_y[mm + nn] = A_nm * (
                                self.H_nm[nn, mm, glob_n] * self.H_nm[nn + 1, mm + 1, glob_n].conj() + self.H_nm[
                            nn, - mm, glob_n] * self.H_nm[nn + 1, - mm - 1, glob_n].conj())
                    f_z[mm + nn] = B_nm * self.H_nm[nn, mm, glob_n] * self.H_nm[nn + 1, mm, glob_n].conj()
                f_nx[nn] = psi * f_x.sum()
                f_ny[nn] = psi * f_y.sum()
                f_nz[nn] = psi * f_z.sum()

            coef = 1 / 8 / math.pi ** 2 / self.wave.rho / self.wave.c ** 2 / self.wave.k ** 2
            force_x = f_nx.sum()
            force_x = coef * force_x.real
            force_y = f_ny.sum()
            force_y = coef * force_y.imag
            force_z = f_nz.sum()
            force_z_acc = f_nz[0:nn - 2].sum()
            accuracy_z[glob_n] = np.abs(force_z.real - force_z_acc.real) / force_z.real

            force_z = -2 * coef * force_z.real
            fz_sum = - np.real(2 * coef * f_nz.copy())

            force[glob_n] = [force_y, force_x, force_z] / self.power * self.wave.c
        return force

    def build_rad_force_old(self, force):
        Fx = np.reshape(force[:, 0], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fy = np.reshape(force[:, 1], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fz = np.reshape(force[:, 2], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])

        x_coord = np.reshape(self.points[:, 0], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        y_coord = np.reshape(self.points[:, 1], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        z_coord = np.reshape(self.points[:, 2], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])

        x_grid = np.array(x_coord[:, 0, 0])
        y_grid = np.array(y_coord[0, :, 0])
        z_grid = np.array(z_coord[0, 0, :])
        if (self.coordinates.nx > 1) and (self.coordinates.nz > 1) and (self.coordinates.ny == 1):
            print('xz')
            fig = make_subplots(rows=1, cols=3, subplot_titles=('Fx', 'Fy', 'Fz'),
                                x_title="z, mm", y_title="x, mm", horizontal_spacing=0.16
                                )
            fig.add_trace(go.Contour(
                z=Fx[:, 0, :],
                x=z_grid * 1e3,
                y=x_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.24),
                zmin=np.min(np.min(Fx[:, 0, :], 0)), zmax=np.max(np.max(Fx[:, 0, :], 0))

            ), 1, 1)
            fig.add_trace(go.Contour(
                z=Fy[:, 0, :],
                x=z_grid * 1e3,
                y=x_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.62),
                zmin=np.min(np.min(Fy[:, 0, :], 0)), zmax=np.max(np.max(Fy[:, 0, :], 0))

            ), 1, 2)
            fig.add_trace(go.Contour(
                z=Fz[:, 0, :],
                x=z_grid * 1e3,
                y=x_grid * 1e3,
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=1),
                zmin=np.min(np.min(Fz[:, 0, :], 0)), zmax=np.max(np.max(Fz[:, 0, :], 0))

            ), 1, 3)

            fig.update_layout(height=400,
                              width=1100)
            fig.show()
            fig2 = None
            fig3 = None

            Z, X = np.meshgrid(z_grid, x_grid)
            Fxz = np.abs(Fx ** 2 + Fz ** 2) ** 0.5
            fig0, ax0 = plt.subplots()
            strm = ax0.streamplot(Z * 1e3, X * 1e3, Fz[:, 0, :], Fx[:, 0, :], density=1, color=Fxz[:, 0, :],
                                  linewidth=1, cmap=plt.cm.viridis)
            fig0.colorbar(strm.lines)
            ax0.set_xlabel('z, mm')
            ax0.set_ylabel('x, mm')
            ax0.set_title('Streamline on xz-plane')
            fig0.set_size_inches(10, 5)

        elif (self.coordinates.ny > 1) and (self.coordinates.nz > 1) and (self.coordinates.nx == 1):
            print('yz')
            Z, Y = np.meshgrid(z_grid, y_grid)
            Fyz = np.abs(Fy ** 2 + Fz ** 2) ** 0.5
            fig0, ax0 = plt.subplots()
            strm = ax0.streamplot(Z * 1e3, Y * 1e3, Fz[0, :, :], Fy[0, :, :], density=1, color=Fyz[0, :, :],
                                  linewidth=1, cmap=plt.cm.viridis)
            fig0.colorbar(strm.lines)
            ax0.set_xlabel('z, mm')
            ax0.set_ylabel('y, mm')
            ax0.set_title('Streamline on yz-plane')
            fig0.set_size_inches(10, 5)

            fig = make_subplots(rows=1, cols=3, subplot_titles=('Fx', 'Fy', 'Fz'),
                                x_title="z, mm", y_title="y, mm", horizontal_spacing=0.16
                                )
            fig.add_trace(go.Contour(
                z=Fx[0, :, :],
                x=z_grid * 1e3,
                y=y_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.24),
                zmin=np.min(np.min(Fx[0, :, :], 0)), zmax=np.max(np.max(Fx[0, :, :], 0))

            ), 1, 1)
            fig.add_trace(go.Contour(
                z=Fy[0, :, :],
                x=z_grid * 1e3,
                y=y_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.62),
                zmin=np.min(np.min(Fy[0, :, :], 0)), zmax=np.max(np.max(Fy[0, :, :], 0))

            ), 1, 2)
            fig.add_trace(go.Contour(
                z=Fz[0, :, :],
                x=z_grid * 1e3,
                y=x_grid * 1e3,
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=1),
                zmin=np.min(np.min(Fz[0, :, :], 0)), zmax=np.max(np.max(Fz[0, :, :], 0))

            ), 1, 3)

            fig.update_layout(height=400,
                              width=1100)
            fig.show()
            fig2 = None
            fig3 = None

        elif (self.coordinates.nx > 1) and (self.coordinates.ny > 1) and (self.coordinates.nz == 1):
            print('xy')
            Y, X = np.meshgrid(y_grid, x_grid)
            Fxy = np.abs(Fx ** 2 + Fz ** 2) ** 0.5
            fig0, ax0 = plt.subplots()
            strm = ax0.streamplot(Y * 1e3, X * 1e3, Fy[:, :, 0], Fx[:, :, 0], density=1, color=Fxy[:, :, 0],
                                  linewidth=1, cmap=plt.cm.viridis)
            fig0.colorbar(strm.lines)
            ax0.set_xlabel('y, mm')
            ax0.set_ylabel('x, mm')
            ax0.set_title('Streamline on xy-plane')
            fig0.set_size_inches(7, 5)

            fig = make_subplots(rows=1, cols=3, subplot_titles=('Fx', 'Fy', 'Fz'),
                                x_title="y, mm", y_title="x, mm", horizontal_spacing=0.16
                                )
            fig.add_trace(go.Contour(
                z=Fx[:, :, 0],
                x=y_grid * 1e3,
                y=x_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.24),
                zmin=np.min(np.min(Fx[:, :, 0], 0)), zmax=np.max(np.max(Fx[:, :, 0], 0))

            ), 1, 1)
            fig.add_trace(go.Contour(
                z=Fy[:, :, 0],
                x=y_grid * 1e3,
                y=x_grid * 1e3,
                line=dict(smoothing=0.85),
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=0.62),
                zmin=np.min(np.min(Fy[:, :, 0], 0)), zmax=np.max(np.max(Fy[:, :, 0], 0))

            ), 1, 2)
            fig.add_trace(go.Contour(
                z=Fz[:, :, 0],
                x=y_grid * 1e3,
                y=x_grid * 1e3,
                colorscale='Viridis',
                contours_coloring='heatmap',
                colorbar=dict(len=1, x=1),
                zmin=np.min(np.min(Fz[:, :, 0], 0)), zmax=np.max(np.max(Fz[:, :, 0], 0))

            ), 1, 3)

            fig.update_layout(height=400,
                              width=1100)
            fig.show()
            fig2 = None
            fig3 = None

        elif (self.coordinates.nx == 1) and (self.coordinates.ny == 1) and (self.coordinates.nz > 1):
            fig, ax = plt.subplots(figsize=(10, 5))
            fig2 = None
            fig3 = None
            fig0 = None
            ax.plot(self.points[:, 2] * 1000, force[:, 0])
            ax.plot(self.points[:, 2] * 1000, force[:, 1])
            ax.plot(self.points[:, 2] * 1000, force[:, 2])
            lgnd = ax.legend(['F_x', 'F_y', 'F_z'], loc='upper right', shadow=True)
            ax.set_xlabel('z, mm')
            ax.set_ylabel('Radiation force, N')
            ax.set_title('Radiation force on z-axis')
        elif (self.coordinates.nx == 1) and (self.coordinates.ny > 1) and (self.coordinates.nz == 1):
            fig, ax = plt.subplots(figsize=(10, 5))
            fig2 = None
            fig3 = None
            fig0 = None
            ax.plot(self.points[:, 1] * 1000, force[:, 0])
            ax.plot(self.points[:, 1] * 1000, force[:, 1])
            ax.plot(self.points[:, 1] * 1000, force[:, 2])
            lgnd = ax.legend(['F_x', 'F_y', 'F_z'], loc='upper right', shadow=True)
            ax.set_xlabel('y, mm')
            ax.set_ylabel('Radiation force, N')
            ax.set_title('Radiation force on y-axis')
        elif (self.coordinates.nx > 1) and (self.coordinates.ny == 1) and (self.coordinates.nz == 1):
            fig, ax = plt.subplots(figsize=(10, 5))
            fig2 = None
            fig3 = None
            fig0 = None
            ax.plot(self.points[:, 0] * 1000, force[:, 0])
            ax.plot(self.points[:, 0] * 1000, force[:, 1])
            ax.plot(self.points[:, 0] * 1000, force[:, 2])
            lgnd = ax.legend(['F_x', 'F_y', 'F_z'], loc='upper right', shadow=True)
            ax.set_xlabel('x, mm')
            ax.set_ylabel('Radiation force, N')
            ax.set_title('Radiation force on x-axis')
        elif (self.coordinates.nx > 1) and (self.coordinates.ny > 1) and (self.coordinates.nz > 1):
            num_steps = len(Fx[0, 0, :])
            fig = go.Figure(go.Contour(z=Fx[:, :, 0],
                                       x=x_grid * 1e3,
                                       y=y_grid * 1e3,
                                       colorscale='Viridis'
                                       )
                            )

            for i in range(1, num_steps):
                fig.add_contour(z=Fx[:, :, i],
                                x=x_grid * 1e3,
                                y=y_grid * 1e3,
                                visible=False,
                                colorscale='Viridis',
                                contours_coloring='heatmap')

            steps = []
            for i in range(num_steps):
                step = dict(
                    label=str(np.around(1e3 * z_grid[i], 2)) + 'mm',
                    method='restyle',
                    args=['visible', [False] * num_steps],
                )
                step['args'][1][i] = True
                steps.append(step)

            sliders = [dict(
                steps=steps,
            )]

            fig.layout.sliders = sliders

            fig.update_layout(title_text="Fx component of radiation force along z-axis",
                              legend_orientation="h",
                              xaxis_title="y, mm",
                              yaxis_title="x, mm",
                              autosize=False,
                              height=600,
                              width=600
                              )
            fig.show()
            fig2 = None
            fig3 = None
            fig0 = None

        return fig, fig2, fig3, fig0

    def build_rad_force(self, force):
        Fx = np.reshape(force[:, 0], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fy = np.reshape(force[:, 1], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fz = np.reshape(force[:, 2], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])

        x_coord = np.reshape(self.points[:, 0], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        y_coord = np.reshape(self.points[:, 1], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        z_coord = np.reshape(self.points[:, 2], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])

        x_grid = np.array(x_coord[:, 0, 0]) * 1e3
        y_grid = np.array(y_coord[0, :, 0]) * 1e3
        z_grid = np.array(z_coord[0, 0, :]) * 1e3

        if (self.coordinates.nx > 1) and (self.coordinates.nz > 1) and (self.coordinates.ny == 1):
            type_field = 'xz'
            type_plot = '2d'
            Fx = Fx[:, 0, :]
            Fy = Fy[:, 0, :]
            Fz = Fz[:, 0, :]
            x_axis = z_grid
            y_axis = x_grid
        elif (self.coordinates.ny > 1) and (self.coordinates.nz > 1) and (self.coordinates.nx == 1):
            type_field = 'yz'
            type_plot = '2d'
            Fx = Fx[0, :, :]
            Fy = Fy[0, :, :]
            Fz = Fz[0, :, :]
            x_axis = z_grid
            y_axis = y_grid
        elif (self.coordinates.nx > 1) and (self.coordinates.ny > 1) and (self.coordinates.nz == 1):
            type_field = 'xy'
            type_plot = '2d'
            Fx = Fx[:, :, 0]
            Fy = Fy[:, :, 0]
            Fz = Fz[:, :, 0]
            x_axis = y_grid
            y_axis = x_grid
        elif (self.coordinates.ny > 1) and (self.coordinates.nz > 1) and (self.coordinates.nx > 1):
            type_field = 'xyz'
            type_plot = '3d'
            x_axis = y_grid
            y_axis = x_grid
        elif (self.coordinates.nx > 1) and (self.coordinates.ny == 1) and (self.coordinates.nz == 1):
            type_field = 'x'
            type_plot = '1d'
            Fx = Fx[:, 0, 0]
            Fy = Fy[:, 0, 0]
            Fz = Fz[:, 0, 0]
            x_axis = x_grid
        elif (self.coordinates.ny > 1) and (self.coordinates.nx == 1) and (self.coordinates.nz == 1):
            type_field = 'y'
            type_plot = '1d'
            Fx = Fx[0, :, 0]
            Fy = Fy[0, :, 0]
            Fz = Fz[0, :, 0]
            x_axis = y_grid
        elif (self.coordinates.nz > 1) and (self.coordinates.ny == 1) and (self.coordinates.nx == 1):
            type_field = 'z'
            type_plot = '1d'
            Fx = Fx[0, 0, :]
            Fy = Fy[0, 0, :]
            Fz = Fz[0, 0, :]
            x_axis = z_grid

        # create figure

        if type_plot == '2d':
            fig = go.Figure()

            button_layer_1_height = 0.9
            button_layer_2_height = 0.8
            button_layer_3_height = 0.7
            button_x_1 = -0.2
            button_x_2 = -0.2
            button_x_3 = -0.2

            fig.add_trace(go.Contour(z=Fx,
                                     x=x_axis,
                                     y=y_axis,
                                     colorscale='Viridis',
                                     # contours_coloring='heatmap',
                                     zmin=np.min(np.min(Fx, 0)),
                                     zmax=np.max(np.max(Fx, 0)),
                                     visible=True,
                                     name='Fx'
                                     ))
            fig.add_trace(go.Contour(z=Fy,
                                     x=x_axis,
                                     y=y_axis,
                                     colorscale='Viridis',
                                     # contours_coloring='heatmap',
                                     zmin=np.min(np.min(Fy, 0)),
                                     zmax=np.max(np.max(Fy, 0)),
                                     visible=False,
                                     name='Fy'
                                     ))
            fig.add_trace(go.Contour(z=Fz,
                                     x=x_axis,
                                     y=y_axis,
                                     colorscale='Viridis',
                                     # contours_coloring='heatmap',
                                     zmin=np.min(np.min(Fz, 0)),
                                     zmax=np.max(np.max(Fz, 0)),
                                     visible=False,
                                     name='Fx'
                                     ))
            fig.add_trace(go.Contour(z=(Fx ** 2 + Fy ** 2 + Fz ** 2) ** 0.5,
                                     x=x_axis,
                                     y=y_axis,
                                     colorscale='Viridis',
                                     # contours_coloring='heatmap',
                                     zmin=np.min(np.min((Fx ** 2 + Fy ** 2 + Fz ** 2) ** 0.5, 0)),
                                     zmax=np.max(np.max((Fx ** 2 + Fy ** 2 + Fz ** 2) ** 0.5, 0)),
                                     visible=False,
                                     name='abs(F)'
                                     ))
            # Update plot sizing
            fig.update_layout(
                height=600,
                title="Fx component of normalised radiation force on  " + type_field + "-plane",
                autosize=True,
                margin_l=200,
                template="plotly_white",
            )
            # Add dropdown
            fig.update_layout(
                updatemenus=[
                    dict(
                        buttons=list([
                            dict(
                                args=["type", "contour"],
                                label="Contour",
                                method="restyle"
                            ),
                            dict(
                                args=["type", "heatmap"],
                                label="Heatmap",
                                method="restyle"
                            )
                        ]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=button_x_1,
                        xanchor="left",
                        y=button_layer_1_height,
                        yanchor="top"
                    ),
                    dict(
                        buttons=list([
                            dict(
                                args=["colorscale", "Viridis"],
                                label="Viridis",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Jet"],
                                label="Jet",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Hot"],
                                label="Hot",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Tropic"],
                                label="Spectral",
                                method="restyle"
                            )
                        ]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=button_x_2,
                        xanchor="left",
                        y=button_layer_2_height,
                        yanchor="top"
                    ),
                    dict(
                        active=0,
                        buttons=list([
                            dict(label="Fx",
                                 method="update",
                                 args=[{"visible": [True, False, False, False]},
                                       {
                                           "title": "Fx component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fy",
                                 method="update",
                                 args=[{"visible": [False, True, False, False]},
                                       {
                                           "title": "Fy component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fz",
                                 method="update",
                                 args=[{"visible": [False, False, True, False]},
                                       {
                                           "title": "Fz component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="abs(F)",
                                 method="update",
                                 args=[{"visible": [False, False, False, True]},
                                       {"title": "Normalised radiation force module on  " + type_field + "-plane"}])]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=button_x_3,
                        xanchor="left",
                        y=button_layer_3_height,
                        yanchor="top"
                    )
                ]
            )
            # #Add annotation
            # fig.update_layout(
            #     annotations=[
            #         dict(text="Displayed variable:", showarrow=False,
            #              xref="paper", x=0, y=button_layer_3_height-0.03,
            #              yref="paper", align="right"),
            #         dict(text="Trace type:", showarrow=False,xref="paper",
            #              x=0.0, y=button_layer_1_height-0.03, yref="paper",
            #              align="right"),
            #         dict(text="Colorscale:", x=0, xref="paper",
            #              y=button_layer_2_height-0.03, yref="paper",
            #              align="right", showarrow=False)
            #     ]
            # )

            fig.update_xaxes(title_text=type_field[1] + "-axis, mm")
            fig.update_yaxes(title_text=type_field[0] + "-axis, mm")

        elif type_plot == '1d':
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=x_axis, y=Fx, name="Fx",
                                     line_shape='spline'))
            fig.add_trace(go.Scatter(x=x_axis, y=Fy, name="Fy",
                                     line_shape='spline'))
            fig.add_trace(go.Scatter(x=x_axis, y=Fz, name="Fz",
                                     line_shape='spline'))

            fig.update_traces(hoverinfo='text+name', mode='lines+markers')
            fig.update_layout(legend=dict(y=0.5, traceorder='reversed', font_size=16))
            fig.update_layout(title='Component of normalised radiation force on ' + type_field[0] + '-axis, mm',
                              xaxis_title=type_field[0] + '-axis, mm',
                              yaxis_title='F, c/W*[N]')

        elif type_plot == '3d':
            num_steps = len(Fx[0, 0, :])
            abs_F = (Fx[:, :, :] ** 2 + Fy[:, :, :] ** 2 + Fz[:, :, :] ** 2) ** 0.5
            button_layer_1_height = 0.8
            button_layer_2_height = 0.7
            button_layer_3_height = 0.6
            button_layer_4_height = 1
            button_x_1 = -0.3
            button_x_2 = -0.3
            button_x_3 = -0.3
            button_x_4 = -0.3

            fig = go.Figure(
                data=[go.Contour(z=Fx[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=True,
                                 name='Fx',
                                 colorscale='Viridis'),
                      go.Contour(z=Fy[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='Fy',
                                 colorscale='Viridis'),
                      go.Contour(z=Fz[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='Fz',
                                 colorscale='Viridis'),
                      go.Contour(z=abs_F[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='abs(F)',
                                 colorscale='Viridis'),
                      go.Heatmap(z=Fx[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='Fx',
                                 colorscale='Viridis'),
                      go.Heatmap(z=Fy[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='Fy',
                                 colorscale='Viridis'),
                      go.Heatmap(z=Fz[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='Fz',
                                 colorscale='Viridis'),
                      go.Heatmap(z=abs_F[:, :, 0],
                                 x=x_grid,
                                 y=y_grid,
                                 visible=False,
                                 name='abs(F)',
                                 colorscale='Viridis')
                      ])
            fig.update_layout(height=600,
                              title="Fx component of normalised radiation force on  " + type_field + "-plane",
                              xaxis_title=type_field[1] + '-axis, mm',
                              yaxis_title=type_field[0] + '-axis, mm')
            frames = []
            for k in range(num_steps):
                frames.append(go.Frame(name=str(k),
                                       data=[go.Contour(z=Fx[:, :, k]),
                                             go.Contour(z=Fy[:, :, k]),
                                             go.Contour(z=Fz[:, :, k]),
                                             go.Contour(z=abs_F[:, :, k]),
                                             go.Heatmap(z=Fx[:, :, k]),
                                             go.Heatmap(z=Fy[:, :, k]),
                                             go.Heatmap(z=Fz[:, :, k]),
                                             go.Heatmap(z=abs_F[:, :, k])]))

            steps = []
            for i in range(num_steps):
                step = dict(
                    label=np.array2string(z_grid[i]),
                    method="animate",
                    args=[[str(i)]]
                )
                steps.append(step)

            sliders = [dict(
                steps=steps,
                currentvalue={"prefix": "z-coordinate, mm: ", "font": {"size": 16}},
            )]

            fig.update_layout(
                margin_l=200,
                updatemenus=[
                    # dict(
                    #     buttons=list([
                    #         dict(
                    #             args=["type", "contour"],
                    #             label="Contour",
                    #             method="restyle"
                    #         ),
                    #         dict(
                    #             args=["type", "heatmap"],
                    #             label="Heatmap",
                    #             method="restyle"
                    #         )
                    #     ]),
                    #     direction="down",
                    #     pad={"r": 10, "t": 10},
                    #     showactive=True,
                    #     x=button_x_1,
                    #     xanchor="left",
                    #     y=button_layer_1_height,
                    #     yanchor="top"
                    # ),
                    dict(
                        buttons=list([
                            dict(
                                args=["colorscale", "Viridis"],
                                label="Viridis",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Jet"],
                                label="Jet",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Hot"],
                                label="Hot",
                                method="restyle"
                            ),
                            dict(
                                args=["colorscale", "Tropic"],
                                label="Spectral",
                                method="restyle"
                            )
                        ]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=button_x_2,
                        xanchor="left",
                        y=button_layer_2_height,
                        yanchor="top"
                    ),
                    dict(
                        active=0,
                        buttons=list([
                            dict(label="Fx Contour",
                                 method="update",
                                 args=[{"visible": [True, False, False, False, False, False, False, False]},
                                       {
                                           "title": "Fx component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fy Contour",
                                 method="update",
                                 args=[{"visible": [False, True, False, False, False, False, False, False]},
                                       {
                                           "title": "Fy component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fz Contour",
                                 method="update",
                                 args=[{"visible": [False, False, True, False, False, False, False, False]},
                                       {
                                           "title": "Fz component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="abs(F) Contour",
                                 method="update",
                                 args=[{"visible": [False, False, False, True, False, False, False, False]},
                                       {"title": "Normalised radiation force module on  " + type_field + "-plane"}]),
                            dict(label="Fx Heatmap",
                                 method="update",
                                 args=[{"visible": [False, False, False, False, True, False, False, False]},
                                       {
                                           "title": "Fx component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fy Heatmap",
                                 method="update",
                                 args=[{"visible": [False, False, False, False, False, True, False, False]},
                                       {
                                           "title": "Fy component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="Fz Heatmap",
                                 method="update",
                                 args=[{"visible": [False, False, False, False, False, False, True, False]},
                                       {
                                           "title": "Fz component of normalised radiation force on  " + type_field + "-plane"}]),
                            dict(label="abs(F) Heatmap",
                                 method="update",
                                 args=[{"visible": [False, False, False, False, False, False, False, True]},
                                       {"title": "Normalised radiation force module on  " + type_field + "-plane"}])
                        ]),
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=button_x_3,
                        xanchor="left",
                        y=button_layer_3_height,
                        yanchor="top"
                    ),
                    dict(showactive=False,
                         type="buttons",
                         buttons=[dict(label="Play",
                                       method="animate",
                                       args=[None, {"fromcurrent": True}]),
                                  dict(label="Pause",
                                       method="animate",
                                       args=[[None], {"frame": {"duration": 0,
                                                                "redraw": False},
                                                      "mode": "immediate",
                                                      "transition": {"duration": 0}}])
                                  ],
                         pad={"r": 10, "t": 10},
                         x=button_x_4,
                         xanchor="left",
                         y=button_layer_4_height,
                         yanchor="top"
                         ),

                ]
            )
            fig.layout.sliders = sliders
            fig.frames = frames

        # fig.show()

        return fig