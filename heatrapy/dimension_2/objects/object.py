"""Contains the class object.

Used to create a thermal object and apply or remove fields.

"""

from ... import mats
import os
import copy
import numpy as np


class Object:
    """Object class.

    This class creates thermal objects to be used in more complex systems.
    It includes two methods to apply and remove fields.

    """

    def __init__(self, amb_temperature, material='Cu', dx=0.01, dy=0.01, dt=0.1,
                 size=(10, 10), file_name=None, boundaries=(0, 0, 0, 0),
                 Q=[], Q0=[], initial_state=False, heat_save=False,
                 materials_path=False):
        """Thermal object initialization.

        amb_temperature: ambient temperature of the whole system
        materials: list of strings of all the used materials present in the
            folder materials
        borders: list of the points where there is a change of material
        materials_order: list of the materials list indexes that defines the
            material properties given by borders
        dx: the space step
        dt: the times step
        file_name: file name where the temperature and heat flux are saved
        boundaries: list of two entries that define the boundary condition
            for temperature. If 0 the boundary condition is insulation
        Q: list of fixed heat source coefficient.
        Q0: list of temperature dependent heat source coefficient.
        initial_state: initial state of the materials. True if applied field
            and False is removed field.
        heat_save: True if saving the heat at the two borders.

        """
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
        cond12 = isinstance(heat_save, bool)
        condition = cond01 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10 and cond11 and cond12
        if not condition:
            raise ValueError

        self.materials = [material]
        self.materials_name = [material]
        self.boundaries = boundaries
        self.amb_temperature = amb_temperature

        if materials_path is False:
            tadi = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'tadi.txt'
            tadd = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'tadd.txt'
            cpa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'cpa.txt'
            cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'cp0.txt'
            k0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'k0.txt'
            ka = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'ka.txt'
            rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'rho0.txt'
            rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'rhoa.txt'
            lheat0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'lheat0.txt'
            lheata = os.path.dirname(os.path.realpath(__file__)) + \
                '/../../database/' + self.materials[0] + '/' + 'lheata.txt'
            self.materials[0] = mats.CalMatPro(
                tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata)
        else:
            tadi = materials_path + self.materials[0] + '/' + 'tadi.txt'
            tadd = materials_path + self.materials[0] + '/' + 'tadd.txt'
            cpa = materials_path + self.materials[0] + '/' + 'cpa.txt'
            cp0 = materials_path + self.materials[0] + '/' + 'cp0.txt'
            k0 = materials_path + self.materials[0] + '/' + 'k0.txt'
            ka = materials_path + self.materials[0] + '/' + 'ka.txt'
            rho0 = materials_path + self.materials[0] + '/' + 'rho0.txt'
            rhoa = materials_path + self.materials[0] + '/' + 'rhoa.txt'
            lheat0 = materials_path + self.materials[0] + '/' + 'lheat0.txt'
            lheata = materials_path + self.materials[0] + '/' + 'lheata.txt'
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
            # if self.temperature[-1] != []:
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
                    self.Cp[-1].append(self.materials[self.materials_index[i][j]].cpa(self.amb_temperature))
                    self.rho[-1].append(self.materials[self.materials_index[i][j]].rhoa(self.amb_temperature))
                    self.k[-1].append(self.materials[self.materials_index[i][j]].ka(self.amb_temperature))
                    self.latent_heat[-1].append(
                        self.materials[self.materials_index[i][j]].lheata()
                    )
                    self.lheat[-1].append([])
                    for lh in self.materials[self.materials_index[i][j]].lheata():
                        if self.temperature[i][j][1] < lh[0] and lh[1] > 0.:
                            self.lheat[-1].append([lh[0], 0.])
                        if self.temperature[i][j][1] > lh[0] and lh[1] > 0.:
                            self.lheat[-1].append([lh[0], lh[1]])
                        if self.temperature[i][j][1] < lh[0] and lh[1] < 0.:
                            self.lheat[-1].append([lh[0], -lh[1]])
                        if self.temperature[i][j][1] > lh[0] and lh[1] < 0.:
                            self.lheat[-1].append([lh[0], 0.])
                else:
                    # print(i,j)
                    # print(self.materials_index)
                    self.Cp[-1].append(self.materials[self.materials_index[i][j]].cp0(self.amb_temperature))
                    self.rho[-1].append(self.materials[self.materials_index[i][j]].rho0(self.amb_temperature))
                    self.k[-1].append(self.materials[self.materials_index[i][j]].k0(self.amb_temperature))
                    self.latent_heat[-1].append(
                        self.materials[self.materials_index[i][j]].lheat0()
                    )
                    self.lheat[-1].append([])
                    for lh in self.materials[self.materials_index[i][j]].lheat0():
                        if self.temperature[i][j][1] < lh[0] and lh[1] > 0.:
                            self.lheat[-1].append([lh[0], 0.])
                        if self.temperature[i][j][1] > lh[0] and lh[1] > 0.:
                            self.lheat[-1].append([lh[0], lh[1]])
                        if self.temperature[i][j][1] < lh[0] and lh[1] < 0.:
                            self.lheat[-1].append([lh[0], -lh[1]])
                        if self.temperature[i][j][1] > lh[0] and lh[1] < 0.:
                            self.lheat[-1].append([lh[0], 0.])

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
            if heat_save:
                line = line + ',Q (J/m)'
            line = line + '\n'
            f = open(self.file_name, 'a')
            f.write(line)
            f.close()

    def activate(self, initial_point, final_point, shape='square'):
        """Activation of the material.

        activates a given space interval of the material,
        between the initial_point and final_point.

        """
        if shape == 'square':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            final_point_x = int(final_point[0])
            final_point_y = int(final_point[1])
            # print(self.Cp)
            # print(self.k)
            # print(self.rho)
            # print(self.lheat)
            # print(self.state)
            # print('\n')
            for i in range(initial_point_x, final_point_x):
                for j in range(initial_point_y, final_point_y):
                    if self.state[i][j] is False:
                        self.temperature[i][j][0] = self.temperature[i][j][0] + \
                            self.materials[self.materials_index[i][j]].tadi(
                                self.temperature[i][j][0])
                        self.rho[i][j] = self.materials[self.materials_index[i][j]].rhoa(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[self.materials_index[i][j]].cpa(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[self.materials_index[i][j]].ka(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[self.materials_index[i][j]].lheata()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                        self.state[i][j] = True
                    else:
                        message = 'point ({:d},{:d}) already activated'.format(i, j)
                        print(message)
            # print(self.temperature)
            # print(self.Cp)
            # print(self.k)
            # print(self.rho)
            # print(self.lheat)
            # print(self.state)
            # print('\n\n\n\n')

        if shape == 'circle':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            radius = int(final_point)
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    length = np.sqrt((i-initial_point_x)**2+(j-initial_point_y)**2)
                    if length <= radius:
                        if self.state[i][j] is False:
                            # print(self.temperature)
                            self.temperature[i][j][0] = self.temperature[i][j][0] + \
                                self.materials[self.materials_index[i][j]].tadi(
                                    self.temperature[i][j][0])
                            # print(self.temperature)
                            self.rho[i][j] = self.materials[self.materials_index[i][j]].rhoa(
                                self.temperature[i][j][0])
                            self.Cp[i][j] = self.materials[self.materials_index[i][j]].cpa(
                                self.temperature[i][j][0])
                            self.k[i][j] = self.materials[self.materials_index[i][j]].ka(
                                self.temperature[i][j][0])
                            self.lheat[i][j] = []
                            valh = self.materials[self.materials_index[i][j]].lheata()
                            self.latent_heat[i][j] = valh
                            for lh in self.latent_heat[i][j]:
                                if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                                if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], lh[1]])
                                if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], -lh[1]])
                                if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                            self.state[i][j] = True
                        else:
                            message = 'point ({:d},{:d}) already activated'.format(i, j)
                            print(message)

    def deactivate(self, initial_point, final_point, shape='square'):
        """Deactivation of the material.

        deactivates a given space interval of the material,
        between the initial_point and final_point.

        """
        if shape == 'square':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            final_point_x = int(final_point[0])
            final_point_y = int(final_point[1])
            for i in range(initial_point_x, final_point_x):
                for j in range(initial_point_y, final_point_y):
                    if self.state[i][j] is True:
                        self.temperature[i][j][0] = self.temperature[i][j][0] - \
                            self.materials[self.materials_index[i][j]].tadd(
                                self.temperature[i][j][0])
                        self.rho[i][j] = self.materials[self.materials_index[i][j]].rho0(
                            self.temperature[i][j][0])
                        self.Cp[i][j] = self.materials[self.materials_index[i][j]].cp0(
                            self.temperature[i][j][0])
                        self.k[i][j] = self.materials[self.materials_index[i][j]].k0(
                            self.temperature[i][j][0])
                        self.lheat[i][j] = []
                        valh = self.materials[self.materials_index[i][j]].lheat0()
                        self.latent_heat[i][j] = valh
                        for lh in self.latent_heat[i][j]:
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                        self.state[i][j] = False
                    else:
                        message = 'point ({:d},{:d}) already deactivated'.format(i, j)
                        print(message)

        if shape == 'circle':
            initial_point_x = int(initial_point[0])
            initial_point_y = int(initial_point[1])
            radius = int(final_point)
            for i in range(self.size[0]):
                for j in range(self.size[1]):
                    length = np.sqrt((i-initial_point_x)**2+(j-initial_point_y)**2)
                    if length <= radius:
                        if self.state[i][j] is True:
                            self.temperature[i][j][0] = self.temperature[i][j][0] - \
                                self.materials[self.materials_index[i][j]].tadd(
                                    self.temperature[i][j][0])
                            self.rho[i][j] = self.materials[self.materials_index[i][j]].rho0(
                                self.temperature[i][j][0])
                            self.Cp[i][j] = self.materials[self.materials_index[i][j]].cp0(
                                self.temperature[i][j][0])
                            self.k[i][j] = self.materials[self.materials_index[i][j]].k0(
                                self.temperature[i][j][0])
                            self.lheat[i][j] = []
                            valh = self.materials[self.materials_index[i][j]].lheat0()
                            self.latent_heat[i][j] = valh
                            for lh in self.latent_heat[i][j]:
                                if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                                if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                    self.lheat[i][j].append([lh[0], lh[1]])
                                if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], -lh[1]])
                                if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                                    self.lheat[i][j].append([lh[0], 0.])
                            self.state[i][j] = False
                        else:
                            message = 'point ({:d},{:d}) already deactivated'.format(i, j)
                            print(message)

    def square(self, material='Gd', initial_point=(3, 3), length=(3, 3),
               state=False, materials_path=False):
        """Activation of the material.

        activates a given space interval of the material,
        between the initial_point and final_point.

        """
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
                tadi = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'tadi.txt'
                tadd = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'tadd.txt'
                cpa = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'cpa.txt'
                cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'cp0.txt'
                k0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'k0.txt'
                ka = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'ka.txt'
                rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'rho0.txt'
                rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'rhoa.txt'
                lheat0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'lheat0.txt'
                lheata = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'lheata.txt'
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))
            else:
                tadi = materials_path + self.materials_name[index] + '/' + 'tadi.txt'
                tadd = materials_path + self.materials_name[index] + '/' + 'tadd.txt'
                cpa = materials_path + self.materials_name[index] + '/' + 'cpa.txt'
                cp0 = materials_path + self.materials_name[index] + '/' + 'cp0.txt'
                k0 = materials_path + self.materials_name[index] + '/' + 'k0.txt'
                ka = materials_path + self.materials_name[index] + '/' + 'ka.txt'
                rho0 = materials_path + self.materials_name[index] + '/' + 'rho0.txt'
                rhoa = materials_path + self.materials_name[index] + '/' + 'rhoa.txt'
                lheat0 = materials_path + self.materials_name[index] + '/' + 'lheat0.txt'
                lheata = materials_path + self.materials_name[index] + '/' + 'lheata.txt'
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
                        if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                            self.lheat[i][j].append([lh[0], 0.])
                        if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                            self.lheat[i][j].append([lh[0], lh[1]])
                        if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                            self.lheat[i][j].append([lh[0], -lh[1]])
                        if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
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
                        if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                            self.lheat[i][j].append([lh[0], 0.])
                        if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                            self.lheat[i][j].append([lh[0], lh[1]])
                        if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                            self.lheat[i][j].append([lh[0], -lh[1]])
                        if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                            self.lheat[i][j].append([lh[0], 0.])

    def circle(self, material='Gd', initial_point=(3, 3), radius=3,
               state=False, materials_path=False):
        """Activation of the material.

        activates a given space interval of the material,
        between the initial_point and final_point.

        """
        initial_point_x = int(initial_point[0])
        initial_point_y = int(initial_point[1])
        if material in self.materials_name:
            index = self.materials_name.index(material)
        else:
            index = len(self.materials)
            self.materials_name.append(material)
            if materials_path is False:
                tadi = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'tadi.txt'
                tadd = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'tadd.txt'
                cpa = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'cpa.txt'
                cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'cp0.txt'
                k0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'k0.txt'
                ka = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'ka.txt'
                rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'rho0.txt'
                rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'rhoa.txt'
                lheat0 = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'lheat0.txt'
                lheata = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../../database/' + self.materials_name[index] + '/' + 'lheata.txt'
                self.materials.append(mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata))
            else:
                tadi = materials_path + self.materials_name[index] + '/' + 'tadi.txt'
                tadd = materials_path + self.materials_name[index] + '/' + 'tadd.txt'
                cpa = materials_path + self.materials_name[index] + '/' + 'cpa.txt'
                cp0 = materials_path + self.materials_name[index] + '/' + 'cp0.txt'
                k0 = materials_path + self.materials_name[index] + '/' + 'k0.txt'
                ka = materials_path + self.materials_name[index] + '/' + 'ka.txt'
                rho0 = materials_path + self.materials_name[index] + '/' + 'rho0.txt'
                rhoa = materials_path + self.materials_name[index] + '/' + 'rhoa.txt'
                lheat0 = materials_path + self.materials_name[index] + '/' + 'lheat0.txt'
                lheata = materials_path + self.materials_name[index] + '/' + 'lheata.txt'
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
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
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
                            if self.temperature[i][j][0] < lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], 0.])
                            if self.temperature[i][j][0] > lh[0] and lh[1] > 0.:
                                self.lheat[i][j].append([lh[0], lh[1]])
                            if self.temperature[i][j][0] < lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], -lh[1]])
                            if self.temperature[i][j][0] > lh[0] and lh[1] < 0.:
                                self.lheat[i][j].append([lh[0], 0.])

    def power_add(self, initial_point, final_point, power, shape='square', power_type='Q'):

        # if self.Q0 == []:
        #     for i in range(self.size[0]):
        #         self.Q0.append([])
        #         for j in range(self.size[1]):
        #             self.Q0[-1].append(0.)
        # if self.Q == []:
        #     for i in range(self.size[0]):
        #         self.Q.append([])
        #         for j in range(self.size[1]):
        #             self.Q[-1].append(0.)

        if power_type == 'Q':
            # if self.Q == []:
            #     for i in range(self.size[0]):
            #         self.Q.append([])
            #         for j in range(self.size[1]):
            #             self.Q[-1].append(0.)

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
                        length = np.sqrt((i-initial_point_x)**2+(j-initial_point_y)**2)
                        if length <= radius:
                            self.Q[i][j] = power

        if power_type == 'Q0':
            # if self.Q0 == []:
            #     for i in range(self.size[0]):
            #         self.Q0.append([])
            #         for j in range(self.size[1]):
            #             self.Q0[-1].append(0.)

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
                        length = np.sqrt((i-initial_point_x)**2+(j-initial_point_y)**2)
                        if length <= radius:
                            self.Q0[i][j] = power

    def power_reset(self, power_type='Q'):
        # if power_type == 'Q':
        self.Q = []
        # if power_type == 'Q0':
        self.Q0 = []
