# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the classes system_objects and single_object.

Used to compute system models

"""

from .. import mats
import os
import copy
from .. import solvers
import object


class system_objects:

    """system_objects class

    This class creates a system of thermal objects, establishes contact between
    them and computes the respective thermal processes.

    """

    def __init__(self, number_objects=2, materials=('Cu', 'Cu'),
                 objects_length=(10, 10), amb_temperature=293, dx=0.01, dt=0.1,
                 file_name='data', initial_state=False,
                 boundaries=((2, 0), (3, 0))):
        """Initializes the object.

        amb_temperature: ambient temperature of the whole system
        materials: list of strings of all the used materials present in the
            folder materials
        number_objects: integer for the number of thermal objects
        objects_length: list of the object lengths (spacial steps)
        dx: the space step
        dt: the times step
        file_name: file name where the temperature and heat flux are saved
        boundaries: tuple of two entries that define the boundary condition
            for tempreture. The first corresponds to the thermal obect while
            the second defines the temperature. If 0 the boundary condition is
            insulation
        initial_state: initial state of the materials. True if applied field
            and False is removed field.

        """

        # check the validity of inputs
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(materials, tuple)
        cond03 = isinstance(number_objects, int)
        cond04 = isinstance(objects_length, tuple)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, unicode)
        cond07 = cond07 or isinstance(file_name, str)
        cond08 = isinstance(boundaries, tuple)
        cond09 = isinstance(initial_state, bool)
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        if not condition:
            raise ValueError

        # initial definitions
        self.objects = []
        for i in range(number_objects):
            if i not in [l[0] for l in boundaries] or (i, 0) in boundaries:
                heat_save = False
            else:
                heat_save = True

            self.objects.append(object.object(amb_temperature,
                                materials=(materials[i],),
                                borders=(1, objects_length[i]+1),
                                materials_order=(0,), dx=dx, dt=dt,
                                file_name=file_name+'_'+str(i)+'.txt',
                                boundaries=(0, 0), Q=[], Q0=[],
                                initial_state=initial_state,
                                heat_save=heat_save))

        self.contacts = set()
        self.boundaries = boundaries
        self.dt = dt
        self.q1 = 0.
        self.q2 = 0.

        for i in boundaries:
            if i[1] != 0:
                for j in range(len(self.objects[i[0]].temperature)):
                    self.objects[i[0]].temperature[j] = [i[1], i[1]]

    def contactFilter(self, object):
        """Filter self.contacts by thermal object id

        object: thermal object id

        """

        filtered = [x for x in
                    self.contacts if (x[0][0] == object or x[1][0] == object)]
        return set(filtered)

        self.contacts.add(contact)

    def contactAdd(self, contact):
        """Add contact to self.contacts

        contact: thermal contact

        """

        self.contacts.add(contact)

    def contactRemove(self, contact):
        """Remove all contacts from an object

        contact: thermal object id

        """

        removing_contact = None

        for i in range(len(self.contacts)):
            if self.contacts[i][0] == contact:
                removing_contact = self.contacts[i]

        self.contacts.remove(removing_contact)

    def compute(self, timeInterval, writeInterval, solver='implicit_k'):
        """Computes the thermal process

        Computes the system for timeInterval, and writes into the file_name
        file every writeInterval time steps. Four different solvers can be
        used: 'explicit_general', 'explicit_k(x)', 'implicit_general',
        and 'implicit_k(x)'.

        """

        # number of time steps for the given timeInterval
        nt = int(timeInterval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for j in range(nt):

            for obj in self.objects:
                obj.Q0 = copy.copy(obj.Q0_ref)

            for contact in self.contacts:
                td1 = self.objects[contact[1][0]].temperature[contact[1][1]][0]
                td2 = self.objects[contact[0][0]].temperature[contact[0][1]][0]
                heat_contact_1 = contact[2] * (td1 - td2)
                heat_contact_2 = contact[2] * (td2 - td1)
                self.objects[contact[0][0]].Q0[contact[0][1]] = heat_contact_1
                self.objects[contact[1][0]].Q0[contact[1][1]] = heat_contact_2

            object_number = -1
            for obj in self.objects:
                object_number = object_number + 1
                obj.time_passed = obj.time_passed + obj.dt

                cond1 = object_number not in [l[0] for l in self.boundaries]
                if cond1 or (object_number, 0) in self.boundaries:

                    # defines the material properties
                    for i in range(1, obj.num_points - 1):
                        if obj.state[i] is True:
                            ind = obj.materials_index[i]
                            obj.rho[i] = obj.materials[ind].rhoa(
                                obj.temperature[i][0])
                            obj.Cp[i] = obj.materials[ind].cpa(
                                obj.temperature[i][0])
                            obj.k[i] = obj.materials[ind].ka(
                                obj.temperature[i][0])
                        if obj.state[i] is False:
                            ind = obj.materials_index[i]
                            obj.rho[i] = obj.materials[ind].rho0(
                                obj.temperature[i][0])
                            obj.Cp[i] = obj.materials[ind].cp0(
                                obj.temperature[i][0])
                            obj.k[i] = obj.materials[ind].k0(
                                obj.temperature[i][0])

                    # SOLVERS

                    # implicit k constant
                    if solver == 'implicit_general':
                        obj.temperature = solvers.implicit_general(obj)

                    # implicit k dependent on x
                    if solver == 'implicit_k(x)':
                        obj.temperature = solvers.implicit_k(obj)

                    # explicit k constant
                    if solver == 'explicit_general':
                        obj.temperature = solvers.explicit_general(obj)

                    # explicit k dependent on x
                    if solver == 'explicit_k(x)':
                        obj.temperature = solvers.explicit_K(obj)

                    # writes the temperature to file_name file ...
                    # if the number of time steps is verified
                    if nw + 1 == writeInterval or j == 0 or j == nt - 1:
                        line = '%f' % obj.time_passed
                        for i in obj.temperature:
                            new_line = ',%f' % i[1]
                            line = line + new_line
                        f = open(obj.file_name, 'a')
                        f.write(line+'\n')
                        f.close()

                else:

                    heat = [p*self.dt*obj.dx for p in obj.Q0 if p is not None]
                    heat = sum(heat)/(len(heat)*obj.dx)

                    if object_number == self.boundaries[0][0]:
                        self.q1 = self.q1 + heat
                        q = self.q1
                    else:
                        self.q2 = self.q2 + heat
                        q = self.q2

                    # writes the temperature to file_name file ...
                    # if the number of time steps is verified
                    if nw + 1 == writeInterval or j == 0 or j == nt - 1:
                        line = '%f' % obj.time_passed
                        for i in obj.temperature:
                            new_line = ',%f' % i[1]
                            line = line + new_line
                        new_line = ',%f' % q
                        line = line + new_line
                        f = open(obj.file_name, 'a')
                        f.write(line+'\n')
                        f.close()

            if nw == writeInterval:
                nw = 0
            else:
                nw = nw + 1


class single_object(object.object):

    """single_object class

    This class solves numerically the heat conduction equation for 1 dimension
    of an active material(s). Three solvers can be used: explicit with
    x-independent k, explicit with x-dependent k, implicit with x-independent
    k, and implicit with x-dependent k. The class has 5 methods: activate for
    the activation of part of the solid, deactivate for the deactivation of
    part of the solid, and compute for solving the equation for a given period
    of time. This class is suited for simulations involving caloric systems
    such as magnetocaloric or electrocaloric systems.

    """

    def __init__(self, amb_temperature, materials=('Cu'), borders=(1, 11),
                 materials_order=(0,), dx=0.01, dt=0.1, file_name='data.txt',
                 boundaries=(0, 0), Q=[], Q0=[], heat_points=(1, -2),
                 initial_state=False, h_left=50000., h_right=50000.):
        """Initializes the object.

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
            for tempreture. If 0 the boundary condition is insulation
        Q: list of 3 entry lists that gives the fixed heat source coeficient.
            The first term is the initial space index where it is applies. The
            second is the final space index where it is applies. The third is
            the value of the coeficient.
        Q0 is a list of 3 entry lists that gives the temperature dependent heat
            source coefficient. The first term is the initial space index where
            it is applies. The second is the final space index where it is
            applies. The third is the value of the coeficient.
        heat_points: list of the space indexes where we want to extract the
            heat flux. Normally, the first term is the heat flux of the hot end
            and the second term is the heat flux of the cold end
        initial_state: initial state of the materials. True if applied field
            and False is removed field.
        h_left: left heat transfer coefficient
        h_right: right heat transfer coefficient

        """

        # check the validity of inputs
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(materials, tuple)
        cond03 = isinstance(borders, tuple)
        cond04 = isinstance(materials_order, tuple)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, unicode)
        cond07 = cond07 or isinstance(file_name, str)
        cond08 = isinstance(boundaries, tuple)
        cond09 = isinstance(heat_points, tuple)
        cond10 = isinstance(initial_state, bool)
        cond11 = isinstance(h_left, int) or isinstance(h_left, float)
        cond12 = isinstance(h_right, int) or isinstance(h_right, float)
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10 and cond11 and cond12
        if not condition:
            raise ValueError

        # initial definitions
        self.heat_points = heat_points
        self.borders = borders
        self.materials = range(len(materials))
        self.boundaries = boundaries
        self.amb_temperature = amb_temperature
        self.h_left = h_left
        self.h_right = h_right

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
        self.temperature.append([amb_temperature, amb_temperature])
        self.rho.append(None)
        self.Cp.append(None)
        self.Q.append(None)
        self.Q0.append(None)
        self.k.append(
            self.materials[self.materials_index[-2]].k0(amb_temperature))

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
        self.state = [initial_state for i in range(borders[-1] + 1)]
        self.time_passed = 0.
        self.heatLeft = 0.
        self.heatRight = 0.

        line = 'time(s)'
        for i in range(len(self.temperature)):
            line = line + ',T[' + str(i) + '] (K)'
        line = line + ',heat[' + str(self.heat_points[0]) + '](W)' + \
            ',heat[' + str(self.heat_points[1]) + '](J)\n'
        f = open(self.file_name, 'a')
        f.write(line)
        f.close()

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

        Computes the system for timeInterval, and writes into the file_name
        file every writeInterval time steps. Four different solvers can be
        used: 'explicit_general', 'explicit_k(x)', 'implicit_general',
        and 'implicit_k(x)'. heat_points is a list that defines the points
        where the heat flux are calculated if modeTemp is True the compute
        method stops when the point modeTempPoint changes numFlag relative to
        the initial value

        """

        # number of time steps for the given timeInterval
        nt = int(timeInterval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # saves the reference temperature for the modeTemp
        initModeTemp = self.temperature[modeTempPoint][0]

        # computes
        for j in range(nt):

            # updates the time_passed
            self.time_passed = self.time_passed + self.dt

            # defines the material properties accoring to the state list
            for i in range(1, self.num_points - 1):
                if self.state[i] is True:
                    self.rho[i] = self.materials[self.materials_index[i]].rhoa(
                        self.temperature[i][0])
                    self.Cp[i] = self.materials[self.materials_index[i]].cpa(
                        self.temperature[i][0])
                    self.k[i] = self.materials[self.materials_index[i]].ka(
                        self.temperature[i][0])
                if self.state[i] is False:
                    self.rho[i] = self.materials[self.materials_index[i]].rho0(
                        self.temperature[i][0])
                    self.Cp[i] = self.materials[self.materials_index[i]].cp0(
                        self.temperature[i][0])
                    self.k[i] = self.materials[self.materials_index[i]].k0(
                        self.temperature[i][0])

            # SOLVERS

            # implicit k constant
            if solver == 'implicit_general':
                self.temperature = solvers.implicit_general(self)

            # implicit k dependent on x
            if solver == 'implicit_k(x)':
                self.temperature = solvers.implicit_k(self)

            # explicit k constant
            if solver == 'explicit_general':
                self.temperature = solvers.explicit_general(self)

            # explicit k dependent on x
            if solver == 'explicit_k(x)':
                self.temperature = solvers.explicit_k(self)

            # calculates the heat flux of the defined ...
            # two points during the time step

            self.heatLeft = (-self.dt * self.h_left *
                             (self.temperature[self.heat_points[0]][1] -
                              self.boundaries[0]) + self.heatLeft)
            self.heatRight = (self.dt * self.h_right *
                              (self.temperature[self.heat_points[1]][1] -
                               self.boundaries[1]) + self.heatRight)

            nw = nw + 1

            # stops the simulation if modeTemp is True ...
            # and the stop condition is verified
            if modeTemp is True:
                value = self.temperature[modeTempPoint][0] - initModeTemp
                if abs(value) > numFlag:
                    nw = 0
                    f = open(self.file_name, 'a')
                    f.write(line)
                    f.close()
                    break

            # writes the temperature to file_name file ...
            # if the number of time steps is verified
            if nw == writeInterval or j == 0 or j == nt - 1:
                line = '%f,' % self.time_passed
                for i in self.temperature:
                    new_line = '%f,' % i[1]
                    line = line + new_line
                new_line = '%f,' % self.heatLeft
                line = line + new_line
                new_line = '%f' % self.heatRight
                line = line + new_line + '\n'
                f = open(self.file_name, 'a')
                f.write(line)
                f.close()

            if nw == writeInterval:
                nw = 0
