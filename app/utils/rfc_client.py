import math
import scipy
import scipy.special as special
import matplotlib.pyplot as plt
import numpy as np
from os.path import join as pjoin
import scipy.io as sio
from scipy.fftpack import fftn, ifftn
import matplotlib.animation as animation

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

    def calculate_force(self):
        force = np.zeros((c.N_points, 3)) * 1.0
        return force

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
        #   self.angle_spectrum = 1/self.spectrum.Nx**2*angle_spectrum_0 * lin_lm *(4*np.pi**2)

        self.angle_spectrum = angle_spectrum_0.conj() * lin_lm * (4 * np.pi ** 2) / self.spectrum.Nx ** 2  #### pi^2???

    #         scipy.io.savemat('ang_spec.mat', {'mydata': self.angle_spectrum})

    #         self.angle_spectrum[19,13]=4*np.pi**2

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
        self.kx_array = self.spectrum.dk / self.spectrum.dx * x_array.copy() - self.spectrum.dk
        self.ky_array = self.spectrum.dk / self.spectrum.dx * y_array.copy() - self.spectrum.dk
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

            force[glob_n] = [force_x, force_y, force_z] / self.power * self.wave.c
        return force  # , scat_p

    def build_rad_force(self, force):
        Fx = np.reshape(force[:, 0], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fy = np.reshape(force[:, 1], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])
        Fz = np.reshape(force[:, 2], [self.coordinates.nx, self.coordinates.ny, self.coordinates.nz])

        if (self.coordinates.nx > 1) and (self.coordinates.nz > 1) and (self.coordinates.ny == 1):
            print('xz')
            fig, ax = plt.subplots()
            extent = [1e3 * min(self.coordinates.z), 1e3 * max(self.coordinates.z), 1e3 * min(self.coordinates.x),
                      1e3 * max(self.coordinates.x)]
            im1 = ax.matshow(Fx[:, 0, :], extent=extent, aspect='auto')
            plt.colorbar(im1)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('x, mm')
            ax.set_title('Fx component of radiation force on xz-plane')
            fig2, ax = plt.subplots()
            im2 = ax.matshow(Fy[:, 0, :], extent=extent, aspect='auto')
            plt.colorbar(im2)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('x, mm')
            ax.set_title('Fy component of radiation force on xz-plane')
            fig3, ax = plt.subplots()
            im3 = ax.matshow(Fz[:, 0, :], extent=extent, aspect='auto')
            plt.colorbar(im3)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('x, mm')
            ax.set_title('Fz component of radiation force on xz-plane')

        elif (self.coordinates.ny > 1) and (self.coordinates.nz > 1) and (self.coordinates.nx == 1):
            print('yz')
            fig, ax = plt.subplots()
            extent = [1e3 * min(self.coordinates.z), 1e3 * max(self.coordinates.z), 1e3 * min(self.coordinates.y),
                      1e3 * max(self.coordinates.y)]
            im1 = ax.matshow(Fx[0, :, :], extent=extent, aspect='auto')
            plt.colorbar(im1)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fx component of radiation force on yz-plane')
            fig2, ax = plt.subplots()
            im2 = ax.matshow(Fy[0, :, :], extent=extent, aspect='auto')
            plt.colorbar(im2)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fy component of radiation force on yz-plane')
            fig3, ax = plt.subplots()
            im3 = ax.matshow(Fz[0, :, :], extent=extent, aspect='auto')
            plt.colorbar(im3)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('z, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fz component of radiation force on yz-plane')
        elif (self.coordinates.nx > 1) and (self.coordinates.ny > 1) and (self.coordinates.nz == 1):
            print('xy')
            fig, ax = plt.subplots()
            extent = [1e3 * min(self.coordinates.x), 1e3 * max(self.coordinates.x), 1e3 * min(self.coordinates.y),
                      1e3 * max(self.coordinates.y)]
            im1 = ax.matshow(Fx[:, :, 0], extent=extent, aspect='auto')
            plt.colorbar(im1)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('x, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fx component of radiation force on xy-plane')
            fig2, ax = plt.subplots()
            im2 = ax.matshow(Fy[:, :, 0], extent=extent, aspect='auto')
            plt.colorbar(im2)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('x, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fy component of radiation force on xy-plane')
            fig3, ax = plt.subplots()
            im3 = ax.matshow(Fz[:, :, 0], extent=extent, aspect='auto')
            plt.colorbar(im3)
            ax.xaxis.tick_bottom()
            ax.set_xlabel('x, mm')
            ax.set_ylabel('y, mm')
            ax.set_title('Fz component of radiation force on xy-plane')

            Y, X = 1e3 * self.spectrum.dx * np.mgrid[-self.coordinates.nx / 2:self.coordinates.nx / 2,
                                            -self.coordinates.ny / 2:self.coordinates.ny / 2]
            Fxy = np.abs(Fx ** 2 + Fy ** 2) ** 0.5
            fig0, ax0 = plt.subplots()
            strm = ax0.streamplot(X, Y, Fx[:, :, 0], Fy[:, :, 0], density=2, color=Fxy[:, :, 0], linewidth=1,
                                  cmap=plt.cm.viridis)
            fig0.colorbar(strm.lines)
            ax0.set_xlabel('x, mm')
            ax0.set_ylabel('y, mm')
            ax0.set_title('Radiation force direction on xy-plane')

        elif (self.coordinates.nx == 1) and (self.coordinates.ny == 1) and (self.coordinates.nz > 1):
            fig, ax = plt.subplots(figsize=(10, 5))
            fig2 = None
            fig3 = None
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
            ax.plot(self.points[:, 0] * 1000, force[:, 0])
            ax.plot(self.points[:, 0] * 1000, force[:, 1])
            ax.plot(self.points[:, 0] * 1000, force[:, 2])
            lgnd = ax.legend(['F_x', 'F_y', 'F_z'], loc='upper right', shadow=True)
            ax.set_xlabel('x, mm')
            ax.set_ylabel('Radiation force, N')
            ax.set_title('Radiation force on x-axis')
        return fig, fig2, fig3