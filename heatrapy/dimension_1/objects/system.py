"""Contains the class system_objects.

Used to compute systems of thermal objects.

"""

import copy
from .. import solvers
from . import Object


class SystemObjects:
    """System_objects class.

    This class creates a system of unidimensional thermal objects, establishes
    contact between them and computes the respective thermal processes.

    """

    def __init__(self, number_objects=2, materials=('Cu', 'Cu'),
                 objects_length=(10, 10), amb_temperature=293, dx=0.01, dt=0.1,
                 file_name=None, initial_state=False,
                 boundaries=((2, 0), (3, 0)), materials_path=False):
        """System object initialization.

        `number_objects` is the integer number of thermal objects. `materials`
        is the list of strings of all the used materials present in
        `material_path`. `amb_temperature` is the ambient temperature of the
        whole system. `object_length` is the list of thermal object lengths
        (spacial steps). `dx` and `dt` are the space and time steps,
        respectively. `file_name` is the file name where the temperature is
        saved. `boundaries` is a list of tuples of length two that define each
        boundary condition for temperature. If 0 the boundary condition is
        insulation. `materials_path` is absolute path of the materials
        database. If false, then the materials database is the standard
        heatrapy database.

        """
        # check the validity of inputs
        materials = tuple(materials)
        objects_length = tuple(objects_length)
        boundaries = tuple(boundaries)
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(materials, tuple)
        cond03 = isinstance(number_objects, int)
        cond04 = isinstance(objects_length, tuple)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, str)
        cond07 = cond07 or (file_name is None)
        cond08 = isinstance(boundaries, tuple)
        cond09 = isinstance(initial_state, bool)
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        if not condition:
            raise ValueError

        # initial definitions
        self.objects = []
        for i in range(number_objects):
            if file_name:
                file_name = file_name + '_' + str(i) + '.txt'

            self.objects.append(Object(amb_temperature,
                                materials=(materials[i],),
                                borders=(1, objects_length[i]+1),
                                materials_order=(0,), dx=dx, dt=dt,
                                file_name=file_name, boundaries=(0, 0),
                                Q=[], Q0=[], initial_state=initial_state,
                                materials_path=materials_path))

        self.contacts = set()
        self.boundaries = boundaries
        self.dt = dt
        self.q1 = 0.
        self.q2 = 0.

        for i in boundaries:
            if i[1] != 0:
                for j in range(len(self.objects[i[0]].temperature)):
                    self.objects[i[0]].temperature[j] = [i[1], i[1]]

    def contact_filter(self, object):
        """Filter self.contacts by thermal object id.

        object: thermal object id

        """
        # check the validity of inputs
        condition = object in range(len(self.objects))
        if not condition:
            raise ValueError

        filtered = [x for x in
                    self.contacts if (x[0][0] == object or x[1][0] == object)]
        return set(filtered)

    def contact_add(self, contact):
        """Add contact to self.contacts.

        The `contact` parameter is a tuple of length 3 (one element for thermal
        object A, one for thermal object B, and one for the heat transfer
        coefficient). Each thermal object element is a tuple of length 2 where
        the first element is the index of the thermal object and the second is
        the spatial point index.

        """
        # check the validity of inputs
        if isinstance(contact, list) or isinstance(contact, tuple):
            if len(contact) == 3:
                condition = True
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        self.contacts.add(contact)

    def contact_remove(self, object_one, object_two):
        """Contact removal.

        Removes all contacts between `object_one` id and `object_two` id.

        """
        # check the validity of inputs
        condition = isinstance(object_one, int)
        condition = condition and isinstance(object_two, int)
        if not condition:
            raise ValueError

        contact_list = list(self.contacts)
        for i in range(len(contact_list)):
            cond_1 = contact_list[i][0][0] == object_one
            cond_1 = cond_1 and contact_list[i][1][0] == object_two
            cond_2 = contact_list[i][0][0] == object_two
            cond_2 = cond_2 and contact_list[i][1][0] == object_one
            if cond_1 or cond_2:
                self.contacts.remove(contact_list[i])

    def change_boundaries(self, object_id, boundaries):
        """Change boundaries.

        Changes the `boundaries` of `object_id`.

        """
        # check the validity of inputs
        condition = isinstance(object_id, int)
        condition = condition and isinstance(boundaries, tuple)
        if condition:
            if len(boundaries) == 2:
                condition = True
            else:
                condition = False
        if not condition:
            raise ValueError

        self.objects[object_id].boundaries = boundaries

    def compute(self, time_interval, write_interval, solver='implicit_k(x)',
                verbose=True):
        """Compute the thermal process.

        Computes the system for `time_interval`, and writes into the
        `file_name` file every `write_interval` time steps. Four different
        solvers can be used: `'explicit_general'`, `'explicit_k(x)'`,
        `'implicit_general'`, and `'implicit_k(x)'`. If `verbose = True`, then
        the progress of the computation is shown.

        """
        # check the validity of inputs
        cond1 = isinstance(time_interval, float)
        cond1 = cond1 or isinstance(time_interval, int)
        cond2 = isinstance(write_interval, int)
        cond3 = isinstance(solver, str)
        cond4 = isinstance(verbose, bool)
        condition = cond1 and cond2 and cond3 and cond4
        if not condition:
            raise ValueError

        # number of time steps for the given timeInterval
        nt = int(time_interval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for j in range(nt):
            for obj in self.objects:
                obj.Q0 = copy.copy(obj.Q0_ref)

            for contact in self.contacts:
                ind1 = int(contact[1][1])
                ind2 = int(contact[0][1])
                td1 = self.objects[contact[1][0]].temperature[ind1][0]
                td2 = self.objects[contact[0][0]].temperature[ind2][0]
                heat_contact_1 = contact[2] * (td1 - td2)
                heat_contact_2 = contact[2] * (td2 - td1)
                self.objects[contact[0][0]].Q0[ind2] = heat_contact_1
                self.objects[contact[1][0]].Q0[ind1] = heat_contact_2

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
                        value = solvers.implicit_general(obj)
                        obj.temperature, obj.lheat = value

                    # implicit k dependent on x
                    if solver == 'implicit_k(x)':
                        obj.temperature, obj.lheat = solvers.implicit_k(obj)

                    # explicit k constant
                    if solver == 'explicit_general':
                        value = solvers.explicit_general(obj)
                        obj.temperature, obj.lheat = value

                    # explicit k dependent on x
                    if solver == 'explicit_k(x)':
                        obj.temperature, obj.lheat = solvers.explicit_k(obj)

                    # writes the temperature to file_name file ...
                    # if the number of time steps is verified
                    if obj.file_name:
                        if nw + 1 == write_interval or j == 0 or j == nt - 1:
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
                    if obj.file_name:
                        if nw + 1 == write_interval or j == 0 or j == nt - 1:
                            line = '%f' % obj.time_passed
                            for i in obj.temperature:
                                new_line = ',%f' % i[1]
                                line = line + new_line
                            new_line = ',%f' % q
                            line = line + new_line
                            f = open(obj.file_name, 'a')
                            f.write(line+'\n')
                            f.close()

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('progress:', int(100*j/nt), '%', end='\r')
            else:
                nw = nw + 1

        if verbose:
            print('Finished simulation')
