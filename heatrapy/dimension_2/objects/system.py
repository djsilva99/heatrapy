"""Contains the classes system_objects and single_object.

Used to compute two-dimensional thermal objects

"""

import copy
from .. import solvers
from . import Object


class SystemObjects:
    """System_objects class.

    This class creates a system of 2D thermal objects, establishes contact
    between them and computes the respective thermal processes.

    """

    def __init__(self, number_objects=2, materials=('Cu', 'Cu'),
                 objects_length=((10, 10), (10, 10)), amb_temperature=293,
                 dx=0.01, dy=0.01, dt=0.1, file_name=None, initial_state=False,
                 boundaries=((0, 0, 0, 0), (0, 0, 0, 0)),
                 materials_path=False):
        """System object initialization.

        `number_objects` is the integer number of thermal objects. `materials`
        is the list of strings of all the used materials present in
        `material_path`. `amb_temperature` is the ambient temperature of the
        whole system. `object_length` is the list of thermal object lengths
        (tuple of spacial steps). `dx` and `dy` are the space steps along the
        x- and y-axis, respectively. dt is the time step. `file_name` is the
        file name where the temperature is saved. `boundaries` is a list of
        four entries that define the boundary condition for temperature (left,
        right, bottom, top). If 0, the boundary condition is insulation.
        `materials_path` is absolute path of the materials database. If false,
        then the materials database is the standard heatrapy database.

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
        cond06 = isinstance(dy, int) or isinstance(dy, float)
        cond07 = isinstance(dt, int) or isinstance(dt, float)
        cond08 = isinstance(file_name, str)
        cond08 = cond08 or (file_name is None)
        cond09 = isinstance(boundaries, tuple)
        cond10 = isinstance(initial_state, bool)
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10
        if not condition:
            raise ValueError

        # initial definitions
        self.objects = []
        file_name_obj = None
        for i in range(number_objects):

            if file_name:
                file_name_obj = file_name + '_' + str(i) + '.txt'

            self.objects.append(Object(amb_temperature, material=materials[i],
                                       dx=dx, dy=dy, dt=dt,
                                       file_name=file_name_obj,
                                       boundaries=boundaries[i], Q=[], Q0=[],
                                       initial_state=initial_state,
                                       materials_path=materials_path))
        self.contacts = set()
        self.boundaries = boundaries
        self.dt = dt
        self.q1 = 0.
        self.q2 = 0.

    def contact_filter(self, object_id):
        """Filter thermal contacts.

        Filter thermal contacts by `object_id`.

        """
        # check the validity of inputs
        condition = isinstance(object_id, int)
        if not condition:
            raise ValueError

        value = object_id
        filtered = [x for x in
                    self.contacts if (x[0][0] == value or x[1][0] == value)]
        return set(filtered)

    def contact_add(self, contact):
        """Add contact to self.contacts.

        `contact` is a thermal contact tuple of length 3, where the first and
        second entries correspond to tuples of the thermal objects and points
        (object_id, (x,y)), and the third entry is the heat transfer
        coefficient.

        """
        # check the validity of inputs
        if isinstance(contact, tuple):
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

        Changes the `boundaries` variable of `object_id`.

        """
        # check the validity of inputs
        condition = isinstance(object_id, int)
        condition = condition and isinstance(boundaries, tuple)
        if condition:
            if len(boundaries) == 4:
                condition = True
            else:
                condition = False
        if not condition:
            raise ValueError

        self.objects[object_id].boundaries = boundaries

    def compute(self, time_interval, write_interval, solver='explicit_k(x)',
                verbose=True):
        """Compute the thermal process.

        Computes the system for `time_interval`, and writes into the
        `file_name` file every `write_interval` time steps. Two different
        solvers can be used: `'explicit_general'` and `'explicit_k(x)'`. If
        verbose = True, then the progress of the computation is shown.

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

        # number of time steps for the given time_interval
        nt = int(time_interval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for k in range(nt):
            for obj in self.objects:
                obj.Q0 = copy.copy(obj.Q0_ref)

            for contact in self.contacts:
                in1_x = int(contact[1][1][0])
                in1_y = int(contact[1][1][1])
                in2_x = int(contact[0][1][0])
                in2_y = int(contact[0][1][1])
                td1 = self.objects[contact[1][0]].temperature[in1_x][in1_y][0]
                td2 = self.objects[contact[0][0]].temperature[in2_x][in2_y][0]
                heat_contact_1 = contact[2] * (td1 - td2)
                heat_contact_2 = contact[2] * (td2 - td1)
                self.objects[contact[0][0]].Q0[in2_x][in2_y] = heat_contact_1
                self.objects[contact[1][0]].Q0[in1_x][in2_y] = heat_contact_2

            for obj in self.objects:
                obj.time_passed = obj.time_passed + obj.dt

                # defines the material properties
                for i in range(obj.size[0]):
                    for j in range(obj.size[1]):
                        if obj.state[i][j] is True:
                            ind = obj.materials_index[i][j]
                            obj.rho[i][j] = obj.materials[ind].rhoa(
                                obj.temperature[i][j][0])
                            obj.Cp[i][j] = obj.materials[ind].cpa(
                                obj.temperature[i][j][0])
                            obj.k[i][j] = obj.materials[ind].ka(
                                obj.temperature[i][j][0])
                        if obj.state[i][j] is False:
                            ind = obj.materials_index[i][j]
                            obj.rho[i][j] = obj.materials[ind].rho0(
                                obj.temperature[i][j][0])
                            obj.Cp[i][j] = obj.materials[ind].cp0(
                                obj.temperature[i][j][0])
                            obj.k[i][j] = obj.materials[ind].k0(
                                obj.temperature[i][j][0])

                # explicit k constant
                if solver == 'explicit_general':
                    temp = []
                    for i in range(obj.size[0]):
                        temp.append([])
                        for j in range(obj.size[1]):
                            temp[-1].append(obj.temperature[i][j][0])
                    obj.temperature, obj.lheat = solvers.explicit_general(obj)
                    temp = []
                    for i in range(obj.size[0]):
                        temp.append([])
                        for j in range(obj.size[1]):
                            temp[-1].append(obj.temperature[i][j][0])

                # explicit k constant
                if solver == 'explicit_k(x)':
                    temp = []
                    for i in range(obj.size[0]):
                        temp.append([])
                        for j in range(obj.size[1]):
                            temp[-1].append(obj.temperature[i][j][0])
                    obj.temperature, obj.lheat = solvers.explicit_k(obj)
                    temp = []
                    for i in range(obj.size[0]):
                        temp.append([])
                        for j in range(obj.size[1]):
                            temp[-1].append(obj.temperature[i][j][0])

                if obj.file_name:
                    if nw + 1 == write_interval or k == 0 or k == nt - 1:
                        line = '%f' % obj.time_passed
                        for i in range(obj.size[0]):
                            for j in range(obj.size[1]):
                                new_line = ',%f' % obj.temperature[i][j][1]
                                line = line + new_line
                        f = open(obj.file_name, 'a')
                        f.write(line+'\n')
                        f.close()

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('progress:', int(100*k/nt), '%', end='\r')
            else:
                nw = nw + 1

        if verbose:
            print('Finished simulation')
