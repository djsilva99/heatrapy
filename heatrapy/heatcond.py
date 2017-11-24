# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class heatcond_activemat_1D.

Used to compute a 1-dimensional model of heat conduction

"""

from numpy import *
import calmatpro
import os


class heatcond_activemat_1D:

    """heatcond_activemat_1D class

    This class solves numerically the heat conduction equation for 1 dimension
    of an active material(s). Three solvers can be used: explicit with
    x-independent k, explicit with x-dependent k, implicit with x-independent
    k, and implicit with x-dependent k. The class has 5 methods: activate for
    the activation of part of the solid, deactivate for the deactivation of
    part of the solid, and compute for solving the equation for a given period
    of time. This class is suited for simulations involving caloric systems
    such as magnetocaloric or electrocaloric systems.

    """

    def __init__(self, ambTemperature, materials=['Cu'], borders=[1, 11],
                 materialsOrder=[0], dx=0.01, dt=0.1, fileName='data.txt',
                 boundaries=[0, 0], Q=[], Q0=[], heatPoints=[1, -2],
                 initialState=False):
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
        Q: list of 3 entry lists that gives the fixed heat source coeficient.
            The first term is the initial space index where it is applies. The
            second is the final space index where it is applies. The third is
            the value of the coeficient.
        Q0 is a list of 3 entry lists that gives the temperature dependent heat
            source coefficient. The first term is the initial space index where
            it is applies. The second is the final space index where it is
            applies. The third is the value of the coeficient.
        heatPoints: list of the space indexes where we want to extract the heat
            flux. Normally, the first term is the heat flux of the hot end and
            the second term is the heat flux of the cold end

        """

        # initial definitions
        self.heatPoints = heatPoints
        self.borders = borders
        self.materials = range(len(materials))
        self.boundaries = boundaries
        self.ambTemperature = ambTemperature

        # loads the data for each material
        for i in range(len(materials)):
            tadi = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'tadi.txt'
            tadd = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'tadd.txt'
            cpa = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'cpa.txt'
            cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'cp0.txt'
            k0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'k0.txt'
            ka = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'ka.txt'
            rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'rho0.txt'
            rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                '/materials/' + materials[i] + '/' + 'rhoa.txt'
            self.materials[i] = calmatpro.calmatpro(
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
                ambTemperature))  # rhoTSHot)
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

        # creates the power sources
        for power in Q:
            for j in range(power[1], power[2]):
                self.Q[j] = power[0]
        for power in Q0:
            for j in range(power[1], power[2]):
                self.Q0[j] = power[0]

        # initializes the state of each point (True if active and False if ...
        # unactive), the time, heat flux from the hot side and heat flux ...
        # from the cold side
        self.state = [initialState for i in range(borders[-1] + 1)]
        self.timePassed = 0.
        self.heatLeft = 0.
        self.heatRight = 0.

        line = 'time(s)'
        for i in range(len(self.temperature)):
            line = line + ',T[' + str(i) + '] (K)'
        line = line + ',heat[' + str(self.heatPoints[0]) + '](W)' + \
            ',heat[' + str(self.heatPoints[1]) + '](J)\n'
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

    def changeHeatPower(self, Q=[], Q0=[]):
        """Heat power source change

        Changes the coeficients for the heat power sources

        """

        for power in Q:
            for j in range(power[1], power[2]):
                self.Q[j] = power[0]
        for power in Q0:
            for j in range(power[1], power[2]):
                self.Q0[j] = power[0]

    def compute(self, timeInterval, writeInterval, solver='explicit_k(x)',
                modeTemp=False, numFlag=0.5, modeTempPoint=1):
        """Computes the thermal process

        Computes the system for timeInterval, and writes into the fileName
        file every writeInterval time steps. Four different solvers can be
        used: 'explicit_general', 'explicit_k(x)', 'implicit_general',
        and 'implicit_k(x)'. heatPoints is a list that defines the points where
        the heat flux are calculated if modeTemp is True the compute method
        stops when the point modeTempPoint changes numFlag relative to the
        initial value

        """

        # number of time steps for the given timeInterval
        nt = int(timeInterval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # saves the reference temperature for the modeTemp
        initModeTemp = self.temperature[modeTempPoint][0]

        # computes
        for j in range(nt):

            # updates the timePassed
            self.timePassed = self.timePassed + self.dt

            # defines the material properties accoring to the state list
            for i in range(1, self.numPoints - 1):
                if self.state[i] is True:
                    self.rho[i] = self.materials[self.materialsIndex[i]].rhoa(
                        self.temperature[i][0])
                    self.Cp[i] = self.materials[self.materialsIndex[i]].cpa(
                        self.temperature[i][0])
                    self.k[i] = self.materials[self.materialsIndex[i]].ka(
                        self.temperature[i][0])
                if self.state[i] is False:
                    self.rho[i] = self.materials[self.materialsIndex[i]].rho0(
                        self.temperature[i][0])
                    self.Cp[i] = self.materials[self.materialsIndex[i]].cp0(
                        self.temperature[i][0])
                    self.k[i] = self.materials[self.materialsIndex[i]].k0(
                        self.temperature[i][0])

            # SOLVERS

            # implicit k constant
            if solver == 'implicit_general':

                # initializes the matrixes for the equation systems
                a = zeros((self.numPoints, self.numPoints))
                b = zeros(self.numPoints)

                # left boundary
                a[0][0] = 1
                if self.boundaries[0] == 0:
                    b[0] = self.temperature[1][0]
                else:
                    b[0] = self.boundaries[0]

                # right boundary
                a[self.numPoints - 1][self.numPoints - 1] = 1
                if self.boundaries[1] == 0:
                    value = self.temperature[self.numPoints - 2][0]
                    b[self.numPoints - 1] = value
                else:
                    b[self.numPoints - 1] = self.boundaries[1]

                # creates the matrixes and solves the equation systems
                for i in range(1, self.numPoints - 1):
                    beta = self.k[i] * self.dt / \
                        (2 * self.rho[i] * self.Cp[i] * self.dx * self.dx)
                    sigma = self.dt / (self.rho[i] * self.Cp[i])

                    a[i][i - 1] = -beta
                    a[i][i] = 1 + 2 * beta - sigma * self.Q[i]
                    a[i][i + 1] = -beta
                    b[i] = (1 - 2 * beta - sigma * self.Q[i]) * \
                        self.temperature[i][0] + \
                        beta * self.temperature[i + 1][0] + \
                        beta * self.temperature[i - 1][0] + 2. * sigma * \
                        (self.Q0[i] - self.Q[i] * self.ambTemperature)

                x = linalg.solve(a, b)

                # updates the temperature list
                for i in range(self.numPoints):
                    self.temperature[i][1] = x[i]
                    self.temperature[i][0] = x[i]

            # implicit k dependent on x
            if solver == 'implicit_k(x)':

                # initializes the matrixes for the equation systems
                a = zeros((self.numPoints, self.numPoints))
                b = zeros(self.numPoints)

                # left boundary
                a[0][0] = 1
                if self.boundaries[0] == 0:
                    b[0] = self.temperature[1][0]
                else:
                    b[0] = self.boundaries[0]

                # right boundary
                a[self.numPoints - 1][self.numPoints - 1] = 1
                if self.boundaries[1] == 0:
                    value = self.temperature[self.numPoints - 2][0]
                    b[self.numPoints - 1] = value
                else:
                    b[self.numPoints - 1] = self.boundaries[1]

                # creates the matrixes and solves the equation systems
                for i in range(1, self.numPoints - 1):
                    gamma = 4. * self.rho[i] * self.Cp[i] * \
                        self.dx * self.dx / self.dt

                    a[i][i - 1] = self.k[i - 1] + self.k[i]
                    a[i][i] = -(gamma + self.k[i + 1] + self.k[i - 1] +
                                2. * self.k[i] - 2 * self.dt * self.dt *
                                self.Q[i])
                    a[i][i + 1] = self.k[i + 1] + self.k[i]
                    b[i] = -(self.k[i + 1] + self.k[i]) * \
                        self.temperature[i + 1][0] + \
                        (-gamma + self.k[i + 1] + self.k[i - 1] + 2. *
                            self.k[i] - 2 * self.dt * self.dt * self.Q[i]) * \
                        self.temperature[i][0] - \
                        (self.k[i - 1] + self.k[i]) * \
                        self.temperature[i - 1][0] - \
                        4. * self.dx * self.dx * \
                        (self.Q0[i] - self.Q[i] * self.ambTemperature)

                x = linalg.solve(a, b)

                # updates the temperature list
                for i in range(self.numPoints):
                    self.temperature[i][1] = x[i]
                    self.temperature[i][0] = x[i]

            # explicit k constant
            if solver == 'explicit_general':

                # computes
                for i in range(1, self.numPoints - 1):

                    alpha = self.dt * \
                        self.k[i] / (self.rho[i] * self.Cp[i] *
                                     self.dx * self.dx)
                    beta = self.dt / (self.rho[i] * self.Cp[i])

                    Tnew = (1 + beta * self.Q[i]) * self.temperature[i][0] + \
                        alpha * (self.temperature[i - 1][0] -
                                 2 * self.temperature[i][0] +
                                 self.temperature[i + 1][0]) + \
                        beta * (self.Q0[i] - self.Q[i] * self.ambTemperature)
                    self.temperature[i][1] = Tnew

                # left boundary for next time step
                if self.boundaries[0] == 0:
                    self.temperature[0][1] = self.temperature[1][0]
                else:
                    self.temperature[0][1] = self.boundaries[0]

                # right boundary for next time step
                if self.boundaries[1] == 0:
                    self.temperature[self.numPoints - 1][1] = \
                        self.temperature[self.numPoints - 2][0]
                else:
                    self.temperature[self.numPoints -
                                     1][1] = self.boundaries[1]

                # updates the temperature list
                for i in range(self.numPoints):
                    self.temperature[i][0] = self.temperature[i][1]

            # explicit k dependent on x
            if solver == 'explicit_k(x)':

                # computes
                for i in range(1, self.numPoints - 1):

                    eta = self.dt / \
                        (2. * self.rho[i] * self.Cp[i] * self.dx * self.dx)
                    beta = self.dt / (self.rho[i] * self.Cp[i])

                    Tnew = (1 + beta * self.Q[i]) * self.temperature[i][0] + \
                        eta * ((self.k[i + 1] + self.k[i]) *
                               self.temperature[i + 1][0] -
                               (self.k[i - 1] + self.k[i + 1] + 2 *
                                self.k[i]) *
                               self.temperature[i][0] +
                               (self.k[i - 1] + self.k[i]) *
                               self.temperature[i - 1][0]) + \
                        beta * (self.Q0[i] - self.Q[i] * self.ambTemperature)

                    self.temperature[i][1] = Tnew

                # left boundary for next time step
                if self.boundaries[0] == 0:
                    self.temperature[0][1] = self.temperature[1][1]
                else:
                    self.temperature[0][1] = self.boundaries[0]

                # right boundary for next time step
                if self.boundaries[1] == 0:
                    self.temperature[self.numPoints - 1][1] = \
                        self.temperature[self.numPoints - 2][1]
                else:
                    self.temperature[self.numPoints -
                                     1][1] = self.boundaries[1]

                # updates the temperature list
                for i in range(self.numPoints):
                    self.temperature[i][0] = self.temperature[i][1]

            # calculates the heat flux of the defined ...
            # two points during the time step
            self.heatLeft = -(self.k[self.heatPoints[0]] *
                              self.dt / self.dx) * \
                (self.temperature[self.heatPoints[0]][1] -
                 self.temperature[self.heatPoints[0] - 1][1]) + \
                self.heatLeft
            self.heatRight = -(self.k[self.heatPoints[1]] *
                               self.dt / self.dx) * \
                (self.temperature[self.heatPoints[1]][1] -
                 self.temperature[self.heatPoints[1] - 1][1]) + \
                self.heatRight

            nw = nw + 1

            # stops the simulation if modeTemp is True ...
            # and the stop condition is verified
            if modeTemp is True:
                value = self.temperature[modeTempPoint][0] - initModeTemp
                if abs(value) > numFlag:
                    nw = 0
                    f = open(self.fileName, 'a')
                    f.write(line)
                    f.close()
                    break

            # writes the temperature to fileName file ...
            # if the number of time steps is verified
            if nw == writeInterval or nw == 1:
                line = '%f,' % self.timePassed
                for i in self.temperature:
                    newLine = '%f,' % i[1]
                    line = line + newLine
                newLine = '%f,' % self.heatLeft
                line = line + newLine
                newLine = '%f' % self.heatRight
                line = line + newLine + '\n'
                f = open(self.fileName, 'a')
                f.write(line)
                f.close()

            if nw == writeInterval:
                nw = 0
