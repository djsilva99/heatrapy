"""Contains the class object.

Used to create unidimensional thermal objects and apply some methods. It does
not include the compute method.

"""

from pathlib import Path
from ... import mats
import os
import copy


class Object:
    """Object class.

    This class creates a unidimensional thermal object. It includes two
    methods to apply and remove fields.

    """

    def __init__(self, amb_temperature, materials=('Cu',), borders=(1, 11),
                 materials_order=(0,), dx=0.01, dt=0.1, file_name=None,
                 boundaries=(0, 0), Q=[], Q0=[], initial_state=False,
                 materials_path=False):
        """Thermal object initialization.

        `amb_temperature` is the ambient temperature of the whole system.
        `materials` is the list of strings of all the used materials present in
        `material_path`. `borders` is a list of the points where there is a
        change of material. `materials_order` is a list of the materials list
        indexes that defines the material properties given by borders. `dx` and
        `dt` are the space and time steps, respectively. `file_name` is the
        file name where the temperature is saved. `boundaries` is a list of two
        entries that define the boundary condition for temperature. If 0 the
        boundary condition is insulation. `initial_state` is the initial state
        of the materials. True if there are an applied field and False if them
        field is absent. `materials_path` is absolute path of the materials
        database. If false, then the materials database is the standard
        heatrapy database. `Q` is a list of fixed heat source coefficient and
        `Q0` is a list of temperature dependent heat source coefficient.

        """
        # check the validity of inputs
        materials = tuple(materials)
        borders = tuple(borders)
        materials_order = tuple(materials_order)
        boundaries = tuple(boundaries)
        Q = list(Q)
        Q0 = list(Q0)
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(materials, tuple)
        cond03 = isinstance(borders, tuple)
        cond04 = isinstance(materials_order, tuple)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, str)
        cond07 = cond07 or (file_name is None)
        cond08 = isinstance(boundaries, tuple)
        cond09 = isinstance(Q, list)
        cond10 = isinstance(Q0, list)
        cond11 = isinstance(initial_state, bool)
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10 and cond11
        if not condition:
            raise ValueError

        # initial definitions
        self.borders = borders
        self.materials = [i for i in range(len(materials))]
        self.boundaries = boundaries
        self.amb_temperature = amb_temperature

        # loads the data for each material
        if materials_path == False:
            for i in range(len(materials)):
                tadi = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + materials[i] + '/' +
                            'tadi.txt')
                tadd = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + materials[i] + '/' +
                            'tadd.txt')
                cpa = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + materials[i] + '/' + 'cpa.txt')
                cp0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                           '/../../database/' + materials[i] + '/' + 'cp0.txt')
                k0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + materials[i] + '/' + 'k0.txt')
                ka = Path(os.path.dirname(os.path.realpath(__file__)) +
                          '/../../database/' + materials[i] + '/' + 'ka.txt')
                rho0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + materials[i] + '/' +
                            'rho0.txt')
                rhoa = Path(os.path.dirname(os.path.realpath(__file__)) +
                            '/../../database/' + materials[i] + '/' +
                            'rhoa.txt')
                lheat0 = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + materials[i] + '/' +
                              'lheat0.txt')
                lheata = Path(os.path.dirname(os.path.realpath(__file__)) +
                              '/../../database/' + materials[i] + '/' +
                              'lheata.txt')
                self.materials[i] = mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata)
        else:
            for i in range(len(materials)):
                tadi = Path(materials_path + materials[i] + '/' + 'tadi.txt')
                tadd = Path(materials_path + materials[i] + '/' + 'tadd.txt')
                cpa = Path(materials_path + materials[i] + '/' + 'cpa.txt')
                cp0 = Path(materials_path + materials[i] + '/' + 'cp0.txt')
                k0 = Path(materials_path + materials[i] + '/' + 'k0.txt')
                ka = Path(materials_path + materials[i] + '/' + 'ka.txt')
                rho0 = Path(materials_path + materials[i] + '/' + 'rho0.txt')
                rhoa = Path(materials_path + materials[i] + '/' + 'rhoa.txt')
                lheat0 = Path(materials_path + materials[i] + '/' +
                              'lheat0.txt')
                lheata = Path(materials_path + materials[i] + '/' +
                                  'lheata.txt')
                self.materials[i] = mats.CalMatPro(
                    tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa, lheat0, lheata)

        # defines which are the properties of each material point
        self.materials_index = [None]
        for i in range(len(borders) - 1):
            for j in range(borders[i], borders[i + 1]):
                self.materials_index.append(materials_order[i])
        self.materials_index.append(None)

        # loads the inputs and creates the lists temperature, Cp, rho and k
        self.num_points = borders[-1] + 1
        self.file_name = file_name
        self.dx = dx
        self.dt = dt
        self.temperature = [[amb_temperature, amb_temperature]]
        self.latent_heat = [None]
        self.lheat = [None]
        self.Cp = [None]
        self.rho = [None]
        self.Q = [None]
        self.Q0 = [None]
        self.k = [self.materials[self.materials_index[1]].k0(amb_temperature)]
        for i in range(1, borders[-1]):
            self.temperature.append([amb_temperature, amb_temperature])
            self.rho.append(self.materials[self.materials_index[i]].rho0(
                amb_temperature))
            self.Cp.append(
                self.materials[self.materials_index[i]].cp0(amb_temperature))
            self.k.append(
                self.materials[self.materials_index[i]].k0(amb_temperature))
            self.Q.append(0.)
            self.Q0.append(0.)
            self.latent_heat.append(
                self.materials[self.materials_index[i]].lheat0()
            )
            self.lheat.append([])
            for lh in self.materials[self.materials_index[i]].lheat0():
                if self.temperature[i][1] < lh[0] and lh[1] > 0.:
                    self.lheat[-1].append([lh[0], 0.])
                if self.temperature[i][1] > lh[0] and lh[1] > 0.:
                    self.lheat[-1].append([lh[0], lh[1]])
                if self.temperature[i][1] < lh[0] and lh[1] < 0.:
                    self.lheat[-1].append([lh[0], -lh[1]])
                if self.temperature[i][1] > lh[0] and lh[1] < 0.:
                    self.lheat[-1].append([lh[0], 0.])
        self.temperature.append([amb_temperature, amb_temperature])
        self.rho.append(None)
        self.Cp.append(None)
        self.Q.append(None)
        self.Q0.append(None)
        self.lheat.append(None)
        self.latent_heat.append(None)
        self.k.append(
            self.materials[self.materials_index[-2]].k0(amb_temperature))

        if Q != []:
            self.Q = Q
        if Q0 != []:
            self.Q0 = Q0

        # initializes the state of each point (True if active and False if ...
        # unactive), the time, heat flux from the hot side and heat flux ...
        # from the cold side
        self.state = [initial_state for i in range(borders[-1] + 1)]
        self.time_passed = 0.

        self.Q_ref = copy.copy(self.Q)
        self.Q0_ref = copy.copy(self.Q0)

        if file_name:
            line = 'time(s)'
            for i in range(len(self.temperature)):
                line = line + ',T[' + str(i) + '] (K)'
            line = line + '\n'
            f = open(self.file_name, 'a')
            f.write(line)
            f.close()

    def activate(self, initial_point, final_point):
        """Activation of the material.

        activates a given space interval of the material, between the
        `initial_point` and `final_point`.

        """
        initial_point = int(initial_point)
        final_point = int(final_point)
        for i in range(initial_point, final_point):

            if self.state[i] is False:
                self.temperature[i][0] = self.temperature[i][0] + \
                    self.materials[self.materials_index[i]].tadi(
                        self.temperature[i][0])
                self.rho[i] = self.materials[self.materials_index[i]].rhoa(
                    self.temperature[i][0])
                self.Cp[i] = self.materials[self.materials_index[i]].cpa(
                    self.temperature[i][0])
                self.k[i] = self.materials[self.materials_index[i]].ka(
                    self.temperature[i][0])
                self.lheat[i] = []
                valh = self.materials[self.materials_index[i]].lheata()
                self.latent_heat[i] = valh
                for lh in self.latent_heat[i]:
                    if self.temperature[i][0] < lh[0] and lh[1] > 0.:
                        self.lheat[i].append([lh[0], 0.])
                    if self.temperature[i][0] > lh[0] and lh[1] > 0.:
                        self.lheat[i].append([lh[0], lh[1]])
                    if self.temperature[i][0] < lh[0] and lh[1] < 0.:
                        self.lheat[i].append([lh[0], -lh[1]])
                    if self.temperature[i][0] > lh[0] and lh[1] < 0.:
                        self.lheat[i].append([lh[0], 0.])
                self.state[i] = True

            else:
                message = 'point %f already activated' % float(i)
                print(message)

    def deactivate(self, initial_point, final_point):
        """Deactivation of the material.

        deactivates a given space interval of the material, between the
        `initial_point` and `final_point`.

        """
        initial_point = int(initial_point)
        final_point = int(final_point)
        for i in range(initial_point, final_point):

            if self.state[i] is True:
                self.temperature[i][0] = self.temperature[i][0] - \
                    self.materials[self.materials_index[i]].tadd(
                        self.temperature[i][0])
                self.rho[i] = self.materials[self.materials_index[i]].rho0(
                    self.temperature[i][0])
                self.Cp[i] = self.materials[self.materials_index[i]].cp0(
                    self.temperature[i][0])
                self.k[i] = self.materials[self.materials_index[i]].k0(
                    self.temperature[i][0])
                self.lheat[i] = []
                valh = self.materials[self.materials_index[i]].lheat0()
                self.latent_heat[i] = valh
                for lh in self.latent_heat[i]:
                    if self.temperature[i][0] < lh[0] and lh[1] > 0.:
                        self.lheat[i].append([lh[0], 0.])
                    if self.temperature[i][0] > lh[0] and lh[1] > 0.:
                        self.lheat[i].append([lh[0], lh[1]])
                    if self.temperature[i][0] < lh[0] and lh[1] < 0.:
                        self.lheat[i].append([lh[0], -lh[1]])
                    if self.temperature[i][0] > lh[0] and lh[1] < 0.:
                        self.lheat[i].append([lh[0], 0.])
                self.state[i] = False

            else:
                message = 'point %f already deactivated' % float(i)
                print(message)
