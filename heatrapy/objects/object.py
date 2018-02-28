# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class object.

Used to create a thermal object and apply or remove fields.

"""

from .. import mats
import os
import copy
from .. import solvers


class object:

    """object class

    This class creates thermal objects to be used in more complex systems.
    It includes 2 methods to apply and remove fields.

    """

    def __init__(self, ambTemperature, materials=['Cu'], borders=[1, 11],
                 materialsOrder=[0], dx=0.01, dt=0.1, fileName='data.txt',
                 boundaries=[0, 0], Q=[], Q0=[], initialState=False,
                 heat_save=False):
        """Initializes the object.

        ambTemperature: ambient temperature of the whole system
        materials: list of strings of all the used materials present in the
            folder materials
        borders: list of the points where there is a change of material
        materialsOrder: list of the materials list indexes that defines the
            material properties given by borders
        dx: the space step
        dt: the times step
        fileName: file name where the temperature and heat flux are saved
        boundaries: list of two entries that define the boundary condition
            for tempreture. If 0 the boundary condition is insulation
        Q: list of fixed heat source coeficient.
        Q0: list of temperature dependent heat source coefficient.
        initialState: initial state of the materials. True if applied field
            and False is removed field.
        heat_save: True if saving the heat at the two borders.

        """

        # initial definitions
        self.borders = borders
        self.materials = range(len(materials))
        self.boundaries = boundaries
        self.ambTemperature = ambTemperature

        # loads the data for each material
        for i in range(len(materials)):
            tadi = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'tadi.txt'
            tadd = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'tadd.txt'
            cpa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'cpa.txt'
            cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'cp0.txt'
            k0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'k0.txt'
            ka = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'ka.txt'
            rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'rho0.txt'
            rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + materials[i] + '/' + 'rhoa.txt'
            self.materials[i] = mats.calmatpro(
                tadi, tadd, cpa, cp0, k0, ka, rho0, rhoa)

        # defines which are the properties of each material point
        self.materialsIndex = [None]
        for i in range(len(borders) - 1):
            for j in range(borders[i], borders[i + 1]):
                self.materialsIndex.append(materialsOrder[i])
        self.materialsIndex.append(None)

        # loads the inputs and creates the lists temperature, Cp, rho and k
        self.numPoints = borders[-1] + 1
        self.fileName = fileName
        self.dx = dx
        self.dt = dt
        self.temperature = [[ambTemperature, ambTemperature]]
        self.Cp = [None]
        self.rho = [None]
        self.Q = [None]
        self.Q0 = [None]
        self.k = [self.materials[self.materialsIndex[1]].k0(ambTemperature)]
        for i in range(1, borders[-1]):
            self.temperature.append([ambTemperature, ambTemperature])
            self.rho.append(self.materials[self.materialsIndex[i]].rho0(
                ambTemperature))
            self.Cp.append(
                self.materials[self.materialsIndex[i]].cp0(ambTemperature))
            self.k.append(
                self.materials[self.materialsIndex[i]].k0(ambTemperature))
            self.Q.append(0.)
            self.Q0.append(0.)
        self.temperature.append([ambTemperature, ambTemperature])
        self.rho.append(None)
        self.Cp.append(None)
        self.Q.append(None)
        self.Q0.append(None)
        self.k.append(
            self.materials[self.materialsIndex[-2]].k0(ambTemperature))

        if Q != []:
            self.Q = Q
        if Q0 != []:
            self.Q0 = Q0

        # initializes the state of each point (True if active and False if ...
        # unactive), the time, heat flux from the hot side and heat flux ...
        # from the cold side
        self.state = [initialState for i in range(borders[-1] + 1)]
        self.timePassed = 0.

        self.Q_ref = copy.copy(self.Q)
        self.Q0_ref = copy.copy(self.Q0)

        line = 'time(s)'
        for i in range(len(self.temperature)):
            line = line + ',T[' + str(i) + '] (K)'
        if heat_save:
            line = line + ',Q (J/m)'
        line = line + '\n'
        f = open(self.fileName, 'a')
        f.write(line)
        f.close()

    def activate(self, initialPoint, finalPoint):
        """Activation of the material

        activates a given space interval of the material,
        between the initialPoint and finalPoint.

        """

        for i in range(initialPoint, finalPoint):

            if self.state[i] is False:
                self.temperature[i][0] = self.temperature[i][0] + \
                    self.materials[self.materialsIndex[i]].tadi(
                        self.temperature[i][0])
                self.rho[i] = self.materials[self.materialsIndex[i]].rhoa(
                    self.temperature[i][0])
                self.Cp[i] = self.materials[self.materialsIndex[i]].cpa(
                    self.temperature[i][0])
                self.k[i] = self.materials[self.materialsIndex[i]].ka(
                    self.temperature[i][0])
                self.state[i] = True

            else:
                message = 'point %f already activated' % float(i)
                print message

    def deactivate(self, initialPoint, finalPoint):
        """Deactivation of the material

        deactivates a given space interval of the material,
        between the initialPoint and finalPoint.

        """

        for i in range(initialPoint, finalPoint):

            if self.state[i] is True:
                self.temperature[i][0] = self.temperature[i][0] - \
                    self.materials[self.materialsIndex[i]].tadd(
                        self.temperature[i][0])
                self.rho[i] = self.materials[self.materialsIndex[i]].rho0(
                    self.temperature[i][0])
                self.Cp[i] = self.materials[self.materialsIndex[i]].cp0(
                    self.temperature[i][0])
                self.k[i] = self.materials[self.materialsIndex[i]].k0(
                    self.temperature[i][0])
                self.state[i] = False

            else:
                message = 'point %f already deactivated' % float(i)
                print message
