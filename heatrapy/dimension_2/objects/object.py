"""Contains the class object.

Used to create two-dimensional thermal objects, and apply some methods. It does
not include the compute method.

"""

from pathlib import Path
from ... import mats
import os
import copy
import numpy as np


class Object:
    """Object class.

    This class creates a two-dimensional thermal object. It includes two
    methods to apply and remove fields.

    """

    def __init__(self, amb_temperature, material='Cu', dx=0.01, dy=0.01,
                 dt=0.1, size=(10, 10), file_name=None,
                 boundaries=(0, 0, 0, 0), Q=[], Q0=[], initial_state=False,
                 materials_path=False):
        """Thermal object initialization.

        `amb_temperature` is the ambient temperature of the whole system.
        `materials` is the background material present in `material_path`.
        `dx`, `dy` are the space steps along the x- and y-axis, respectively.
        `dt` is the time step. `file_name` is the file name where the
        temperature is saved. `boundaries` is a list of four entries that
        define the boundary condition for temperature (left, right, bottom,
        top). If 0 the boundary condition is insulation. `initial_state` is the
        initial state of the materials. True if there are an applied field and
        False if them field is absent. `materials_path` is absolute path of the
        materials database. If false, then the materials database is the
        standard heatrapy database. `Q` is a list of fixed heat source
        coefficient and `Q0` is a list of temperature dependent heat source
        coefficient.

        """
        # check the validity of inputs
        boundaries = tuple(boundaries)
        Q = list(Q)
        Q0 = list(Q0)
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, str)
        cond07 = cond07 or (file_name is None)
        cond08 = isinstance(boundaries, tuple)
        cond09 = isinstance(Q, list)
        cond10 = isinstance(Q0, list)
        cond11 = isinstance(initial_state, bool)
        condition = cond01 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10 and cond11
        if not condition:
            raise ValueError

        self.materials = [material]
        self.materials_name = [material]
        self.boundaries = boundaries
        self.amb_temperature = amb_temperature

        if materials_path is False:
            tadi = Path(os.path.dirname(os.path.realpath(__file__)) +
                        '/../../database/' + self.materials[0] + '/' +
                        'tadi.txt')
            tadd = Path(os.path.dirname(os.path.realpath(__file__)) +
                        '/../../database/' + self.materials[0] + '/' +
                        'tadd.txt')
            cpa = Path(os.path.dirname(os.path.realpath(__file__)) +
                       '/../../database/' + self.materials[0] + '/' +
                       'cpa.txt')
            cp0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                       '/../../database/' + self.materials[0] + '/' +
                       'cp0.txt')
            k0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                      '/../../database/' + self.materials[0] + '/' + 'k0.txt')
            ka = Path(os.path.dirname(os.path.realpath(__file__)) +
                      '/../../database/' + self.materials[0] + '/' + 'ka.txt')
            rho0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                        '/../../database/' + self.materials[0] + '/' +
                        'rho0.txt')
            rhoa = Path(os.path.dirname(os.path.realpath(__file__)) +
                        '/../../database/' + self.materials[0] + '/' +
                        'rhoa.txt')
            lheat0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + self.materials[0] + '/' +
                          'lheat0.txt')
            lheata = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + self.materials[0] + '/' +
                          'lheata.txt')
            self.materials[0] = mats.CalMatPro(
                tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata)
        else:
            tadi = Path(materials_path + self.materials[0] + '/' + 'tadi.txt')
            tadd = Path(materials_path + self.materials[0] + '/' + 'tadd.txt')
            cpa = Path(materials_path + self.materials[0] + '/' + 'cpa.txt')
            cp0 = Path(materials_path + self.materials[0] + '/' + 'cp0.txt')
            k0 = Path(materials_path + self.materials[0] + '/' + 'k0.txt')
            ka = Path(materials_path + self.materials[0] + '/' + 'ka.txt')
            rho0 = Path(materials_path + self.materials[0] + '/' + 'rho0.txt')
            rhoa = Path(materials_path + self.materials[0] + '/' + 'rhoa.txt')
            lheat0 = Path(materials_path + self.materials[0] + '/' +
                          'lheat0.txt')
            lheata = Path(materials_path + self.materials[0] + '/' +
                          'lheata.txt')
            self.materials[0] = mats.CalMatPro(
                tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata)

        self.materials_index = [None]
        self.size = size
        self.file_name = file_name
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.temperature = []
        self.latent_heat = []
        self.lheat = []
        self.Cp = []
        self.rho = []
        self.Q = []
        self.Q0 = []
        self.k = []
        self.materials_index = []
        self.state = []
        for i in range(self.size[0]):
            self.state.append([])
            self.materials_index.append([])
            self.temperature.append([])
            self.latent_heat.append([])
            self.lheat.append([])
            self.Cp.append([])
            self.rho.append([])
            self.Q.append([0. for i in range(self.size[1])])
            self.Q0.append([0. for i in range(self.size[1])])
            self.k.append([])
            for j in range(self.size[1]):
                self.materials_index[-1].append(0)
                self.temperature[-1].append([amb_temperature, amb_temperature])
                self.state[-1].append(initial_state)
                if initial_state:
                    value = self.materials[self.materials_index[i][j]]
                    self.Cp[-1].append(value.cpa(self.amb_temperature))
                    value = self.materials[self.materials_index[i][j]]
                    self.rho[-1].append(value.rhoa(self.amb_temperature))
                    value = self.materials[self.materials_index[i][j]]
                    self.k[-1].append(value.ka(self.amb_temperature))
                    self.latent_heat[-1].append(
                        self.materials[self.materials_index[i][j]].lheata()
                    )
                    self.lheat[-1].append([])
                    value = self.materials[self.materials_index[i][j]]
                    for lh in value.lheata():
                        if self.temperature[i][j][1] < lh[0] and lh[1] > 0.:
                            self.lheat[-1][-1].append([lh[0], 0.])
                        if self.temperature[i][j][1] > lh[0] and lh[1] > 0.:
                            self.lheat[-1][-1].append([lh[0], lh[1]])
                        if self.temperature[i][j][1] < lh[0] and lh[1] < 0.:
                            self.lheat[-1][-1].append([lh[0], -lh[1]])
                        if self.temperature[i][j][1] > lh[0] and lh[1] < 0.:
                            self.lheat[-1][-1].append([lh[0], 0.])
                else:
                    value = self.materials[self.materials_index[i][j]]
                    self.Cp[-1].append(value.cp0(self.amb_temperature))
                    value = self.materials[self.materials_index[i][j]]
                    self.rho[-1].append(value.rho0(self.amb_temperature))
                    value = self.materials[self.materials_index[i][j]]
                    self.k[-1].append(value.k0(self.amb_temperature))
                    self.latent_heat[-1].append(
                        self.materials[self.materials_index[i][j]].lheat0()
                    )
                    self.lheat[-1].append([])
                    value = self.materials[self.materials_index[i][j]]
                    for lh in value.lheat0():
                        if self.temperature[i][j][1] < lh[0] and lh[1] > 0.:
                            self.lheat[-1][-1].append([lh[0], 0.])
                        if self.temperature[i][j][1] > lh[0] and lh[1] > 0.:
                            self.lheat[-1][-1].append([lh[0], lh[1]])
                        if self.temperature[i][j][1] < lh[0] and lh[1] < 0.:
                            self.lheat[-1][-1].append([lh[0], -lh[1]])
                        if self.temperature[i][j][1] > lh[0] and lh[1] < 0.:
                            self.lheat[-1][-1].append([lh[0], 0.])

        if Q != []:
            self.Q = Q
        if Q0 != []:
            self.Q0 = Q0

        self.time_passed = 0.

        self.Q_ref = copy.copy(self.Q)
        self.Q0_ref = copy.copy(self.Q0)

        if file_name:
            line = 'time(s)'
            for i in range(size[0]):
                for j in range(size[1]):
                    line = line + ',T[' + str(i) + '][' + str(j) + '] (K)'
            line = line + '\n'
            f = open(self.file_name, 'a')
            f.write(line)
            f.close()

    def activate(self, initial_point, final_point, shape='square'):
        """Activation of the material.

        Activates a given piece of material. If `shape` is `'square'`, then the
        `initial_point` is the tuple (x,y) of the bottom left point and the
        `final_point` is the tuple (x,y) of the top right point. If the shape
        is `'circle'`, the `initial_point` is the tuple (x,y) of the center of
        the circle and `final_point` is its radius.

        """
        # check the validity of inputs
        if isinstance(shape, str):
            if shape == 'square':
                value = isinstance(initial_point, tuple)
                if value and isinstance(final_point, tuple):
                    condition = len(initial_point) == 2
                    condition = condition and len(final_point) == 2
                else:
                    condition = False
            elif shape == 'circle':
                value = isinstance(final_point, int)
                value = value or isinstance(final_point, float)
                value = value and isinstance(initial_point, tuple)
                if value:
                    condition = len(initial_point) == 2
                else:
                    condition = False
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        if shape == 'square':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            final_point_x = int(final_point[0])
            final_point_y = int(final_point[1])
            for i in range(initial_point_x, final_point_x):
                for j in range(initial_point_y, final_point_y):
                    if self.state[i][j] is False:
                        value = self.temperature[i][j][0]
                        self.temperature[i][j][0] = value + \
                            self.materials[self.materials_index[i][j]].tadi(
                                self.temperature[i][j][0])
                        value = self.materials_index[i][j]
                        self.rho[i][j] = self.materials[value].rhoa(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[value].cpa(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[value].ka(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[value].lheata()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            cond = self.temperature[i][j][0] < lh[0]
                            if cond and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            cond = self.temperature[i][j][0] > lh[0]
                            if cond and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            cond = self.temperature[i][j][0] < lh[0]
                            if cond and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            cond = self.temperature[i][j][0] > lh[0]
                            if cond and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                        self.state[i][j] = True
                    else:
                        message = 'point ({:d},{:d})'.format(i, j)
                        message = message + ' already activated'
                        print(message)

        if shape == 'circle':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            radius = int(final_point)
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    value = (i-initial_point_x)**2+(j-initial_point_y)**2
                    length = np.sqrt(value)
                    if length <= radius:
                        if self.state[i][j] is False:
                            value = self.temperature[i][j][0]
                            index_value = self.materials_index[i][j]
                            self.temperature[i][j][0] = value + \
                                self.materials[index_value].tadi(
                                    self.temperature[i][j][0])
                            self.rho[i][j] = self.materials[index_value].rhoa(
                                self.temperature[i][j][0])
                            self.Cp[i][j] = self.materials[index_value].cpa(
                                self.temperature[i][j][0])
                            self.k[i][j] = self.materials[index_value].ka(
                                self.temperature[i][j][0])
                            self.lheat[i][j] = []
                            valh = self.materials[index_value].lheata()
                            self.latent_heat[i][j] = valh
                            for lh in self.latent_heat[i][j]:
                                value = self.temperature[i][j][0] < lh[0]
                                if value and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                                value = self.temperature[i][j][0] > lh[0]
                                if value and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], lh[1]])
                                value = self.temperature[i][j][0] < lh[0]
                                if value and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], -lh[1]])
                                value = self.temperature[i][j][0] > lh[0]
                                if value and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                            self.state[i][j] = True
                        else:
                            message = 'point ({:d},{:d})'.format(i, j)
                            message = message + ' already activated'
                            print(message)

    def deactivate(self, initial_point, final_point, shape='square'):
        """Deactivation of the material.

        Deactivates a given piece of material. If `shape` is `'square'`, then
        the `initial_point` is the tuple (x,y) of the bottom left point and the
        `final_point` is the tuple (x,y) of the top right point. If the shape
        is `'circle'`, the `initial_point` is the tuple (x,y) of the center of
        the circle and `final_point` is its radius.

        """
        # check the validity of inputs
        if isinstance(shape, str):
            if shape == 'square':
                value = isinstance(initial_point, tuple)
                if value and isinstance(final_point, tuple):
                    condition = len(initial_point) == 2
                    condition = condition and len(final_point) == 2
                else:
                    condition = False
            elif shape == 'circle':
                value = isinstance(final_point, int)
                value = value or isinstance(final_point, float)
                value = value and isinstance(initial_point, tuple)
                if value:
                    condition = len(initial_point) == 2
                else:
                    condition = False
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        if shape == 'square':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            final_point_x = int(final_point[0])
            final_point_y = int(final_point[1])
            for i in range(initial_point_x, final_point_x):
                for j in range(initial_point_y, final_point_y):
                    if self.state[i][j] is True:
                        value = self.temperature[i][j][0]
                        self.temperature[i][j][0] = value - \
                            self.materials[self.materials_index[i][j]].tadd(
                                self.temperature[i][j][0])
                        value = self.materials_index[i][j]
                        self.rho[i][j] = self.materials[value].rho0(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[value].cp0(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[value].k0(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[value].lheat0()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            cond = self.temperature[i][j][0] < lh[0]
                            if cond and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            cond = self.temperature[i][j][0] > lh[0]
                            if cond and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            cond = self.temperature[i][j][0] < lh[0]
                            if cond and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            cond = self.temperature[i][j][0] > lh[0]
                            if cond and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                        self.state[i][j] = False
                    else:
                        message = 'point ({:d},{:d})'.format(i, j)
                        message = message + ' already deactivated'
                        print(message)

        if shape == 'circle':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            radius = int(final_point)
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    value = (i-initial_point_x)**2+(j-initial_point_y)**2
                    length = np.sqrt(value)
                    if length <= radius:
                        if self.state[i][j] is False:
                            value = self.temperature[i][j][0]
                            index_value = self.materials_index[i][j]
                            self.temperature[i][j][0] = value - \
                                self.materials[index_value].tadd(
                                    self.temperature[i][j][0])
                            self.rho[i][j] = self.materials[index_value].rho0(
                                self.temperature[i][j][0])
                            self.Cp[i][j] = self.materials[index_value].cp0(
                                self.temperature[i][j][0])
                            self.k[i][j] = self.materials[index_value].k0(
                                self.temperature[i][j][0])
                            self.lheat[i][j] = []
                            valh = self.materials[index_value].lheat0()
                            self.latent_heat[i][j] = valh
                            for lh in self.latent_heat[i][j]:
                                value = self.temperature[i][j][0] < lh[0]
                                if value and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                                value = self.temperature[i][j][0] > lh[0]
                                if value and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], lh[1]])
                                value = self.temperature[i][j][0] < lh[0]
                                if value and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], -lh[1]])
                                value = self.temperature[i][j][0] > lh[0]
                                if value and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                            self.state[i][j] = False
                        else:
                            message = 'point ({:d},{:d})'.format(i, j)
                            message = message + ' already deactivated'
                            print(message)

    def square(self, material='Gd', initial_point=(3, 3), length=(3, 3),
               state=False, materials_path=False):
        """Material adding with rectangle shape.

        Adds a new material with a rectangle shape, where `initial_point` is
        the bottom left (x,y) tuple, `length` is the length along the two axis,
        state is the initial state of the material and `materials_path` is the
        absolute path of the materials database.

        """
        # check the validity of inputs
        value = isinstance(initial_point, tuple)
        if value and isinstance(length, tuple):
            cond1 = len(initial_point) == 2
            cond1 = cond1 and len(length) == 2
        else:
            cond1 = False
        cond2 = isinstance(material, str)
        cond3 = isinstance(state, bool)
        cond4 = isinstance(materials_path, str) or materials_path is False
        if not cond1 and cond2 and cond3 and cond4:
            raise ValueError

        initial_point_x = int(initial_point[0])
        initial_point_y = int(initial_point[1])
        final_point_x = initial_point_x + int(length[0])
        final_point_y = initial_point_y + int(length[1])
        if material in self.materials_name:
            index = self.materials_name.index(material)
        else:
            index = len(self.materials)
            self.materials_name.append(material)
            if materials_path is False:
                value = self.materials_name[index]
                tadi = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'tadi.txt')
                tadd = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'tadd.txt')
                cpa = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + value + '/' + 'cpa.txt')
                cp0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + value + '/' + 'cp0.txt')
                k0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + value + '/' + 'k0.txt')
                ka = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + value + '/' + 'ka.txt')
                rho0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'rho0.txt')
                rhoa = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'rhoa.txt')
                lheat0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + value + '/' + 'lheat0.txt')
                lheata = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + value + '/' + 'lheata.txt')
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))
            else:
                value = self.materials_name[index]
                tadi = Path(materials_path + value + '/' + 'tadi.txt')
                tadd = Path(materials_path + value + '/' + 'tadd.txt')
                cpa = Path(materials_path + value + '/' + 'cpa.txt')
                cp0 = Path(materials_path + value + '/' + 'cp0.txt')
                k0 = Path(materials_path + value + '/' + 'k0.txt')
                ka = Path(materials_path + value + '/' + 'ka.txt')
                rho0 = Path(materials_path + value + '/' + 'rho0.txt')
                rhoa = Path(materials_path + value + '/' + 'rhoa.txt')
                lheat0 = Path(materials_path + value + '/' + 'lheat0.txt')
                lheata = Path(materials_path + value + '/' + 'lheata.txt')
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))

        for i in range(initial_point_x, final_point_x):
            for j in range(initial_point_y, final_point_y):
                if state is False:
                    self.state[i][j] = False
                    self.materials_index[i][j] = index
                    self.rho[i][j] = self.materials[index].rho0(
                        self.temperature[i][j][0])
                    self.Cp[i][j] = self.materials[index].cp0(
                        self.temperature[i][j][0])
                    self.k[i][j] = self.materials[index].k0(
                        self.temperature[i][j][0])
                    self.lheat[i][j] = []
                    valh = self.materials[index].lheat0()
                    self.latent_heat[i][j] = valh
                    for lh in self.latent_heat[i][j]:
                        if self.temperature[i][j][0] < lh[0] and lh[1] > 0:
                            self.lheat[i][j].append([lh[0], 0.])
                        if self.temperature[i][j][0] > lh[0] and lh[1] > 0:
                            self.lheat[i][j].append([lh[0], lh[1]])
                        if self.temperature[i][j][0] < lh[0] and lh[1] < 0:
                            self.lheat[i][j].append([lh[0], -lh[1]])
                        if self.temperature[i][j][0] > lh[0] and lh[1] < 0:
                            self.lheat[i][j].append([lh[0], 0.])

                else:
                    self.state[i][j] = True
                    self.materials_index[i][j] = index
                    self.rho[i][j] = self.materials[index].rhoa(
                        self.temperature[i][j][0])
                    self.Cp[i][j] = self.materials[index].cpa(
                        self.temperature[i][j][0])
                    self.k[i][j] = self.materials[index].ka(
                        self.temperature[i][j][0])
                    self.lheat[i][j] = []
                    valh = self.materials[index].lheata()
                    self.latent_heat[i][j] = valh
                    for lh in self.latent_heat[i][j]:
                        if self.temperature[i][j][0] < lh[0] and lh[1] > 0:
                            self.lheat[i][j].append([lh[0], 0.])
                        if self.temperature[i][j][0] > lh[0] and lh[1] > 0:
                            self.lheat[i][j].append([lh[0], lh[1]])
                        if self.temperature[i][j][0] < lh[0] and lh[1] < 0:
                            self.lheat[i][j].append([lh[0], -lh[1]])
                        if self.temperature[i][j][0] > lh[0] and lh[1] < 0:
                            self.lheat[i][j].append([lh[0], 0.])

    def circle(self, material='Gd', initial_point=(3, 3), radius=3,
               state=False, materials_path=False):
        """Material adding with circle shape.

        Adds a new material with a circle shape, where `initial_point` is the
        (x,y) tuple of the center of the circle, `radius` is the radius of the
        circle, state is the initial state of the material and `materials_path`
        is the absolute path of the materials database.

        """
        # check the validity of inputs
        if isinstance(initial_point, tuple):
            cond1 = len(initial_point) == 2
        else:
            cond1 = False
        cond2 = isinstance(radius, int) or isinstance(radius, float)
        cond3 = isinstance(material, str)
        cond4 = isinstance(state, bool)
        cond5 = isinstance(materials_path, str) or materials_path is False
        if not cond1 and cond2 and cond3 and cond4 and cond5:
            raise ValueError

        initial_point_x = int(initial_point[0])
        initial_point_y = int(initial_point[1])
        if material in self.materials_name:
            index = self.materials_name.index(material)
        else:
            index = len(self.materials)
            self.materials_name.append(material)
            if materials_path is False:
                value = self.materials_name[index]
                tadi = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'tadi.txt')
                tadd = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'tadd.txt')
                cpa = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + value + '/' + 'cpa.txt')
                cp0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + value + '/' + 'cp0.txt')
                k0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + value + '/' + 'k0.txt')
                ka = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + value + '/' + 'ka.txt')
                rho0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'rho0.txt')
                rhoa = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + value + '/' + 'rhoa.txt')
                lheat0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + value + '/' + 'lheat0.txt')
                lheata = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + value + '/' + 'lheata.txt')
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))
            else:
                value = self.materials_name[index]
                tadi = Path(materials_path + value + '/' + 'tadi.txt')
                tadd = Path(materials_path + value + '/' + 'tadd.txt')
                cpa = Path(materials_path + value + '/' + 'cpa.txt')
                cp0 = Path(materials_path + value + '/' + 'cp0.txt')
                k0 = Path(materials_path + value + '/' + 'k0.txt')
                ka = Path(materials_path + value + '/' + 'ka.txt')
                rho0 = Path(materials_path + value + '/' + 'rho0.txt')
                rhoa = Path(materials_path + value + '/' + 'rhoa.txt')
                lheat0 = Path(materials_path + value + '/' + 'lheat0.txt')
                lheata = Path(materials_path + value + '/' + 'lheata.txt')
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))

        for i in range(self.size[0]):
            for j in range(self.size[1]):
                length = np.sqrt((i-initial_point_x)**2+(j-initial_point_y)**2)
                if length <= radius:
                    if state is False:
                        self.state[i][j] = False
                        self.materials_index[i][j] = index
                        self.rho[i][j] = self.materials[index].rho0(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[index].cp0(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[index].k0(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[index].lheat0()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0:
                                self.lheat[i][j].append([lh[0], 0.])

                    else:
                        self.state[i][j] = True
                        self.materials_index[i][j] = index
                        self.rho[i][j] = self.materials[index].rhoa(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[index].cpa(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[index].ka(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[index].lheata()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0:
                                self.lheat[i][j].append([lh[0], 0.])

    def power_add(self, initial_point, final_point, power, shape='square',
                  power_type='Q'):
        """Power adding.

        Adds a power matrix to the thermal object. If `shape` is `'square'`,
        then the `initial_point` is the tuple (x,y) of the bottom left point
        and the `final_point` is the tuple (x,y) of the top right point. If the
        `shape` is `'circle'`, the `initial_point` is the tuple (x,y) of the
        center of the circle and `final_point` is its radius. `power` is the
        value of the power to add, and `power_type` is the type of power to be
        introduced, which has the value `'Q'` if it is temperature dependent
        and `'Q0'` if it is temperature independent.

        """
        # check the validity of inputs
        if isinstance(shape, str):
            if shape == 'square':
                value = isinstance(initial_point, tuple)
                if value and isinstance(final_point, tuple):
                    cond1 = len(initial_point) == 2
                    cond1 = cond1 and len(final_point) == 2
                else:
                    cond1 = False
            elif shape == 'circle':
                value = isinstance(final_point, int)
                value = value or isinstance(final_point, float)
                value = value and isinstance(initial_point, tuple)
                if value:
                    cond1 = len(initial_point) == 2
                else:
                    cond1 = False
            else:
                cond1 = False
        else:
            cond1 = False
        cond2 = isinstance(power, int) or isinstance(power, float)
        if isinstance(power_type, str):
            if power_type == 'Q' or power_type == 'Q0':
                cond3 = True
            else:
                cond3 = False
        else:
            cond3 = False
        if not cond1 and cond2 and cond3:
            raise ValueError

        if power_type == 'Q':
            if shape == 'square':
                initial_point_x = int(initial_point[0])
                initial_point_y = int(initial_point[1])
                final_point_x = int(final_point[0])
                final_point_y = int(final_point[1])
                for i in range(initial_point_x, final_point_x):
                    for j in range(initial_point_y, final_point_y):
                        self.Q[i][j] = power
            if shape == 'circle':
                initial_point_x = int(initial_point[0])
                initial_point_y = int(initial_point[1])
                radius = int(final_point)
                for i in range(self.size[0]):
                    for j in range(self.size[1]):
                        value = (i-initial_point_x)**2+(j-initial_point_y)**2
                        length = np.sqrt(value)
                        if length <= radius:
                            self.Q[i][j] = power

        if power_type == 'Q0':
            if shape == 'square':
                initial_point_x = int(initial_point[0])
                initial_point_y = int(initial_point[1])
                final_point_x = int(final_point[0])
                final_point_y = int(final_point[1])
                for i in range(initial_point_x, final_point_x):
                    for j in range(initial_point_y, final_point_y):
                        self.Q0[i][j] = power
            if shape == 'circle':
                initial_point_x = int(initial_point[0])
                initial_point_y = int(initial_point[1])
                radius = int(final_point)
                for i in range(self.size[0]):
                    for j in range(self.size[1]):
                        value = (i-initial_point_x)**2+(j-initial_point_y)**2
                        length = np.sqrt(value)
                        if length <= radius:
                            self.Q0[i][j] = power

    def power_reset(self, power_type='Q'):
        """Power reset.

        Resets the power matrix with `power_type` `'Q'` or `'Q0'`, which
        corresponds to the power temperature dependent and temperature
        independent, respectively.

        """
        # check the validity of inputs
        if isinstance(power_type, str):
            if power_type == 'Q' or power_type == 'Q0':
                condition = True
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        self.Q = []
        self.Q0 = []
