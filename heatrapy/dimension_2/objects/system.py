"""Contains the classes system_objects and single_object.

Used to compute system models

"""

from ... import mats
import os
import copy
from .. import solvers
from . import Object
import numpy as np
# import time
import matplotlib.pyplot as plt
import matplotlib

class SystemObjects:
    """System_objects class.

    This class creates a system of thermal objects, establishes contact between
    them and computes the respective thermal processes.

    """

    def __init__(self, number_objects=2, materials=('Cu', 'Cu'),
                 objects_length=((10, 10), (10, 10)), amb_temperature=293, dx=0.01, dt=0.1,
                 file_name=None, initial_state=False,
                 boundaries=((0, 0, 0, 0), (0, 0, 0, 0)), materials_path=False):
        """System object initialization.

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
            if i not in [l[0] for l in boundaries] or (i, 0) in boundaries:
                heat_save = False
            else:
                heat_save = True

            if file_name:
                file_name = file_name + '_' + str(i) + '.txt'

            self.objects.append(Object(amb_temperature,
                                material=materials[i],
                                dx=dx, dt=dt,
                                file_name=file_name, boundaries=boundaries[i],
                                Q=[], Q0=[], initial_state=initial_state,
                                heat_save=heat_save,
                                materials_path=materials_path))

        self.contacts = set()
        self.boundaries = boundaries
        self.dt = dt
        self.q1 = 0.
        self.q2 = 0.

        # for i in boundaries:
        #     if i[1] != 0:
        #         for j in range(len(self.objects[i[0]].temperature)):
        #             self.objects[i[0]].temperature[j] = [i[1], i[1]]

    def contact_filter(self, object_id):
        """Filter self.contacts by thermal object id.

        object: thermal object id

        """
        filtered = [x for x in
                    self.contacts if (x[0][0] == object_id or x[1][0] == object_id)]
        return set(filtered)

    def contact_add(self, contact):
        """Add contact to self.contacts.

        contact: thermal contact

        """
        self.contacts.add(contact)

    def contact_remove(self, object_one, object_two):
        """Remove all contacts from an object.

        contact: thermal object id

        """
        contact_list = list(self.contacts)
        for i in range(len(contact_list)):
            # print(contact_list[i][0][0])
            # print(contact_list[i][1][0])
            # print()
            cond_1 = contact_list[i][0][0] == object_one and contact_list[i][1][0] == object_two
            cond_2 = contact_list[i][0][0] == object_two and contact_list[i][1][0] == object_one
            # print(cond_1)
            # print(cond_2)
            if cond_1 or cond_2:
                # print(contact_list[i])
                # removing_contact = contact_list[i]
                self.contacts.remove(contact_list[i])

    def compute(self, time_interval, write_interval, solver='implicit_k(x)',
                verbose=True):
        """Compute the thermal process.

        Computes the system for timeInterval, and writes into the file_name
        file every write_interval time steps. Four different solvers can be
        used: 'explicit_general', 'explicit_k(x)', 'implicit_general',
        and 'implicit_k(x)'.

        """
        # number of time steps for the given timeInterval
        nt = int(time_interval / self.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for k in range(nt):
            for obj in self.objects:
                obj.Q0 = copy.copy(obj.Q0_ref)

            for contact in self.contacts:
                ind1_x = int(contact[1][1][0])
                ind1_y = int(contact[1][1][1])
                ind2_x = int(contact[0][1][0])
                ind2_y = int(contact[0][1][1])
                td1 = self.objects[contact[1][0]].temperature[ind1_x][ind1_y][0]
                td2 = self.objects[contact[0][0]].temperature[ind2_x][ind2_y][0]
                heat_contact_1 = contact[2] * (td1 - td2)
                heat_contact_2 = contact[2] * (td2 - td1)
                self.objects[contact[0][0]].Q0[ind2_x][ind2_y] = heat_contact_1
                self.objects[contact[1][0]].Q0[ind1_x][ind2_y] = heat_contact_2

            # object_number = -1
            for obj in self.objects:
                # object_number = object_number + 1
                obj.time_passed = obj.time_passed + obj.dt

                # cond1 = object_number not in [l[0] for l in self.boundaries]
                # if cond1 or (object_number, 0) in self.boundaries:

                # defines the material properties
                for i in range(self.size[0]):
                    for j in range(self.size[1]):
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
                # if obj.file_name:
                #     if nw + 1 == write_interval or j == 0 or j == nt - 1:
                #         line = '%f' % obj.time_passed
                #         for i in obj.temperature:
                #             new_line = ',%f' % i[1]
                #             line = line + new_line
                #         f = open(obj.file_name, 'a')
                #         f.write(line+'\n')
                #         f.close()
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

                # else:

                #     heat = [p*self.dt*obj.dx for p in obj.Q0 if p is not None]
                #     heat = sum(heat)/(len(heat)*obj.dx)

                #     if object_number == self.boundaries[0][0]:
                #         self.q1 = self.q1 + heat
                #         q = self.q1
                #     else:
                #         self.q2 = self.q2 + heat
                #         q = self.q2

                #     # writes the temperature to file_name file ...
                #     # if the number of time steps is verified
                #     if obj.file_name:
                #         if nw + 1 == write_interval or j == 0 or j == nt - 1:
                #             line = '%f' % obj.time_passed
                #             for i in obj.temperature:
                #                 new_line = ',%f' % i[1]
                #                 line = line + new_line
                #             new_line = ',%f' % q
                #             line = line + new_line
                #             f = open(obj.file_name, 'a')
                #             f.write(line+'\n')
                #             f.close()

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('progress:', int(100*k/nt), '%', end='\r')
            else:
                nw = nw + 1

        if verbose:
            print('Finished simulation')


class SingleObject(Object):
    """Single_object class.

    This class solves numerically the heat conduction equation for 1 dimension
    of an active material(s). Three solvers can be used: explicit with
    x-independent k, explicit with x-dependent k, implicit with x-independent
    k, and implicit with x-dependent k. The class has 5 methods: activate for
    the activation of part of the solid, deactivate for the deactivation of
    part of the solid, and compute for solving the equation for a given period
    of time. This class is suited for simulations involving caloric systems
    such as magnetocaloric or electrocaloric systems.

    """

    def __init__(self, amb_temperature, material='Cu', dx=0.01, dy=0.01,
                 dt=0.1, size=(10, 10), file_name=None,
                 boundaries=(0, 0, 0, 0), Q=[], Q0=[], initial_state=False,
                 materials_path=False, draw=['temperature','state','materials','Q','Q0'], draw_scale=None):
        """Object initialization.

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
        boundaries = tuple(boundaries)
        Q = list(Q)
        Q0 = list(Q0)
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(material, str)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dt, int) or isinstance(dt, float)
        cond07 = isinstance(file_name, str)
        cond07 = cond07 or (file_name is None)
        cond08 = isinstance(boundaries, tuple)
        cond10 = isinstance(initial_state, bool)
        cond13 = isinstance(Q, list)
        cond14 = isinstance(Q0, list)
        condition = cond01 and cond02 and cond05
        condition = condition and cond06 and cond07 and cond08
        condition = condition and cond10 and cond13
        condition = condition and cond14
        if not condition:
            raise ValueError

        self.materials_path = materials_path

        self.time_passed = 0.
        self.size = size
        self.dt=dt
        self.dx=dx
        self.dy=dy
        self.object = Object(amb_temperature,
                             material=material,
                             dx=dx, dy=dy, dt=dt, size=size,
                             file_name=file_name, boundaries=boundaries,
                             Q=[], Q0=[], initial_state=initial_state,
                             materials_path=materials_path)

        self.draw = draw
        self.draw_scale = draw_scale
        if self.draw:
            cmap = matplotlib.cm.RdBu
            name = 'my_cmap_r'
            reverse = []
            k = []
            for key in cmap._segmentdata:
                k.append(key)
                channel = cmap._segmentdata[key]
                data = []
                for t in channel:
                    data.append((1-t[0],t[2],t[1]))
                reverse.append(sorted(data))
            LinearL = dict(zip(k,reverse))
            my_cmap_r = matplotlib.colors.LinearSegmentedColormap(name, LinearL)
            cmap_r = my_cmap_r  # reverse_colourmap(cmap)
            for drawing in self.draw:
                if drawing == 'temperature':
                    self.figure = plt.figure()
                    self.ax = self.figure.add_subplot(111)
                    temp = []
                    for i in range(self.object.size[0]):
                        temp.append([])
                        for j in range(self.object.size[1]):
                            temp[-1].append(self.object.temperature[i][j][0])
                    if not self.draw_scale:
                        vmax = max(max(temp, key=max))
                        vmin = min(min(temp, key=min))
                        temp = np.array(temp)
                        self.im = self.ax.imshow(temp, vmax=vmax, vmin=vmin, cmap=cmap_r,extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower',interpolation='hamming')
                    else:
                        temp = np.array(temp)
                        self.im = self.ax.imshow(temp, vmax=self.draw_scale[0], vmin=self.draw_scale[1], cmap=cmap_r,extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower',interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax.figure.colorbar(self.im, ax=self.ax, **cbar_kw)
                    self.ax.set_title('Temperature (K)')
                    self.ax.set_xlabel('x axis (m)')
                    self.ax.set_ylabel('y axis (m)')
                    plt.show(block=False)
                if drawing == 'state':
                    self.figure_state = plt.figure()
                    self.ax_state = self.figure_state.add_subplot(111)
                    vmax = 1
                    vmin = 0
                    state = np.array(self.object.state)
                    self.im_state = self.ax_state.imshow(state, vmax=1.5, vmin=-0.5, cmap=plt.get_cmap("gray", 2), extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower')
                    cbar_kw = {}
                    cbarlabel = ""  # "state"
                    qrates = np.array(['active', 'inactive'])
                    norm = matplotlib.colors.BoundaryNorm(np.linspace(0, 1, 2), 1)
                    fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
                    cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
                    cbar_state = self.ax_state.figure.colorbar(self.im_state, ax=self.ax_state, **cbar_kw)
                    cbar_state.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                    self.ax_state.set_title('State')
                    self.ax_state.set_xlabel('x axis (m)')
                    self.ax_state.set_ylabel('y axis (m)')
                    plt.show(block=False)
                if drawing == 'materials' and len(self.object.materials_name) > 1:
                    self.figure_materials = plt.figure()
                    self.ax_materials = self.figure_materials.add_subplot(111)
                    vmax = len(self.object.materials)-1
                    vmin = 0
                    material_id = np.array(self.object.materials_index)
                    self.im_materials = self.ax_materials.imshow(material_id, vmax=vmax, vmin=vmin, cmap=plt.get_cmap("gray", vmax), extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower')
                    cbar_kw = {}
                    cbarlabel = ""
                    qrates = np.array(self.object.materials_name)
                    norm = matplotlib.colors.BoundaryNorm(np.linspace(0, 1, len(self.object.materials)), len(self.object.materials)-1)
                    fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
                    cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
                    cbar_materials = self.ax_materials.figure.colorbar(self.im_materials, ax=self.ax_materials, **cbar_kw)
                    cbar_materials.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                    self.ax_materials.set_title('Materials')
                    self.ax_materials.set_xlabel('x axis (m)')
                    self.ax_materials.set_ylabel('y axis (m)')
                    plt.show(block=False)
                if drawing == 'Q':
                    self.figure_Q = plt.figure()
                    self.ax_Q = self.figure_Q.add_subplot(111)
                    temp = np.array(self.object.Q)
                    vmax = 1
                    vmin = 0
                    self.im_Q = self.ax_Q.imshow(temp, vmax=vmax, vmin=vmin, cmap='inferno',extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower',interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax_Q.figure.colorbar(self.im_Q, ax=self.ax_Q, **cbar_kw)
                    self.ax_Q.set_title('Q (W/m²)')
                    self.ax_Q.set_xlabel('x axis (m)')
                    self.ax_Q.set_ylabel('y axis (m)')
                    plt.show(block=False)
                if drawing == 'Q0':
                    self.figure_Q0 = plt.figure()
                    self.ax_Q0 = self.figure_Q0.add_subplot(111)
                    temp = np.array(self.object.Q0)
                    vmax = 1
                    vmin = 0
                    self.im_Q0 = self.ax_Q0.imshow(temp, vmax=vmax, vmin=vmin, cmap='inferno',extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower',interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax_Q0.figure.colorbar(self.im_Q0, ax=self.ax_Q0, **cbar_kw)
                    self.ax_Q0.set_title('Q0 (W/m²)')
                    self.ax_Q0.set_xlabel('x axis (m)')
                    self.ax_Q0.set_ylabel('y axis (m)')
                    plt.show(block=False)

    def show_figure(self, figure_type):
        if figure_type == 'temperature':
            self.ax = self.figure.add_subplot(111)
            temp = []
            for i in range(self.object.size[0]):
                temp.append([])
                for j in range(self.object.size[1]):
                    temp[-1].append(self.object.temperature[i][j][0])
            if not self.draw_scale:
                vmax = max(max(temp, key=max))
                vmin = min(min(temp, key=min))
                temp = np.array(temp)
                self.im = self.ax.imshow(temp, vmax=vmax, vmin=vmin, cmap='jet',extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower')
            else:
                temp = np.array(temp)
                self.im = self.ax.imshow(temp, vmax=self.draw_scale[0], vmin=self.draw_scale[1], cmap='jet',extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower')
            cbar_kw = {}
            cbarlabel = "temperature (K)"
            cbar = self.ax.figure.colorbar(self.im, ax=self.ax, **cbar_kw)
            cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
            self.ax.set_title('Temperature')
            self.ax.set_xlabel('x axis (m)')
            self.ax.set_ylabel('y axis (m)')
            plt.show(block=False)
        if figure_type == 'state':
            self.figure_state = plt.figure()
            self.ax_state = self.figure_state.add_subplot(111)
            vmax = 1
            vmin = 0
            state = np.array(self.object.state)
            self.im_state = self.ax_state.imshow(state, vmax=1.5, vmin=-0.5, cmap=plt.get_cmap("gray", 2), extent =[0, size[0]*dx, 0, size[1]*dy],origin ='lower')
            cbar_kw = {}
            cbarlabel = ""  # "state"
            qrates = np.array(['active', 'inactive'])
            norm = matplotlib.colors.BoundaryNorm(np.linspace(0, 1, 2), 1)
            fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
            cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
            cbar_state = self.ax_state.figure.colorbar(self.im_state, ax=self.ax_state, **cbar_kw)
            cbar_state.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
            self.ax_state.set_title('State')
            self.ax_state.set_xlabel('x axis (m)')
            self.ax_state.set_ylabel('y axis (m)')
            plt.show(block=False)
        if figure_type == 'materials' and len(self.object.materials_name) > 1:
            self.figure_materials = plt.figure()
            self.ax_materials = self.figure_materials.add_subplot(111)
            vmax = len(self.object.materials)-1
            vmin = 0
            material_id = np.array(self.object.materials_index)
            self.im_materials = self.ax_materials.imshow(material_id, vmax=vmax, vmin=vmin, cmap=plt.get_cmap("PiYG", vmax+1), extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower')
            cbar_kw = {}
            cbarlabel = ""
            materials_name_list = copy.deepcopy(self.object.materials_name)
            materials_name_list.reverse()
            qrates = np.array(materials_name_list)
            norm = matplotlib.colors.BoundaryNorm(np.linspace(0, len(self.object.materials)-1, len(self.object.materials)), len(self.object.materials)-1)
            fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
            cbar_kw = dict(ticks=np.arange(0, len(self.object.materials)+1), format=fmt)
            cbar_materials = self.ax_materials.figure.colorbar(self.im_materials, ax=self.ax_materials, **cbar_kw)
            cbar_materials.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
            self.ax_materials.set_title('Materials')
            self.ax_materials.set_xlabel('x axis (m)')
            self.ax_materials.set_ylabel('y axis (m)')
            plt.show(block=False)
        if figure_type == 'Q':
            self.figure_Q = plt.figure()
            self.ax_Q = self.figure_Q.add_subplot(111)
            temp = np.array(self.object.Q)
            vmax = 1
            vmin = 0
            self.im_Q = self.ax_Q.imshow(temp, vmax=vmax, vmin=vmin, cmap='inferno',extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower',interpolation='hamming')
            cbar_kw = {}
            cbar = self.ax_Q.figure.colorbar(self.im_Q, ax=self.ax_Q, **cbar_kw)
            self.ax_Q.set_title('Q (W/m²)')
            self.ax_Q.set_xlabel('x axis (m)')
            self.ax_Q.set_ylabel('y axis (m)')
            plt.show(block=False)
        if figure_type == 'Q0':
            self.figure_Q0 = plt.figure()
            self.ax_Q0 = self.figure_Q0.add_subplot(111)
            temp = np.array(self.object.Q0)
            vmax = 1
            vmin = 0
            self.im_Q0 = self.ax_Q0.imshow(temp, vmax=vmax, vmin=vmin, cmap='inferno',extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower',interpolation='hamming')
            cbar_kw = {}
            cbar = self.ax_Q0.figure.colorbar(self.im_Q0, ax=self.ax_Q0, **cbar_kw)
            self.ax_Q0.set_title('Q0 (W/m²)')
            self.ax_Q0.set_xlabel('x axis (m)')
            self.ax_Q0.set_ylabel('y axis (m)')
            plt.show(block=False)

    def activate(self, initial_point, final_point, shape='square'):
        self.object.activate(initial_point=initial_point, final_point=final_point, shape=shape)
        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    temp = []
                    for i in range(self.object.size[0]):
                        temp.append([])
                        for j in range(self.object.size[1]):
                            temp[-1].append(self.object.temperature[i][j][0])
                    # print(temp)
                    self.im.set_array(temp)
                    if not self.draw_scale:
                        vmax = max(max(temp, key=max))
                        vmin = min(min(temp, key=min))
                        self.im.set_clim(vmin=vmin)
                        self.im.set_clim(vmax=vmax)
                    self.figure.canvas.draw()
                if drawing == 'state':
                    self.im_state.set_array(self.object.state)
                    self.figure_state.canvas.draw()

    def deactivate(self, initial_point, final_point, shape='square'):
        self.object.deactivate(initial_point=initial_point, final_point=final_point, shape=shape)
        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    temp = []
                    for i in range(self.object.size[0]):
                        temp.append([])
                        for j in range(self.object.size[1]):
                            temp[-1].append(self.object.temperature[i][j][0])
                    # print(temp)
                    self.im.set_array(temp)
                    if not self.draw_scale:
                        vmax = max(max(temp, key=max))
                        vmin = min(min(temp, key=min))
                        self.im.set_clim(vmin=vmin)
                        self.im.set_clim(vmax=vmax)
                    self.figure.canvas.draw()
                if drawing == 'state':
                    self.im_state.set_array(self.object.state)
                    self.figure_state.canvas.draw()

    def change_boundaries(self, boundaries):
        self.object.boundaries = boundaries

    def change_material(self, shape, material, initial_point, length, state):
        if shape == 'square':
            self.object.square(material=material,
                               initial_point=initial_point, length=length,
                               state=state, materials_path=self.materials_path)
        if shape == 'circle':
            self.object.circle(material=material,
                               initial_point=initial_point, radius=length,
                               state=state, materials_path=self.materials_path)
        try:
            plt.close(self.figure_materials)
            # self.figure_materials.close()
        except:
            pass

        if self.draw:
            if 'materials' in self.draw:
                self.figure_materials = plt.figure()
                self.ax_materials = self.figure_materials.add_subplot(111)
                vmax = len(self.object.materials)-1
                vmin = 0
                material_id = np.array(self.object.materials_index)
                self.im_materials = self.ax_materials.imshow(material_id, vmax=vmax, vmin=vmin, cmap=plt.get_cmap("PiYG", vmax+1), extent =[0, self.size[0]*self.dx, 0, self.size[1]*self.dy],origin ='lower')
                cbar_kw = {}
                cbarlabel = ""
                materials_name_list = copy.deepcopy(self.object.materials_name)
                materials_name_list.reverse()
                qrates = np.array(materials_name_list)
                norm = matplotlib.colors.BoundaryNorm(np.linspace(0, len(self.object.materials)-1, len(self.object.materials)), len(self.object.materials)-1)
                fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: qrates[::-1][norm(x)])
                cbar_kw = dict(ticks=np.arange(0, len(self.object.materials)+1), format=fmt)
                cbar_materials = self.ax_materials.figure.colorbar(self.im_materials, ax=self.ax_materials, **cbar_kw)
                cbar_materials.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                self.ax_materials.set_title('Materials')
                self.ax_materials.set_xlabel('x axis (m)')
                self.ax_materials.set_ylabel('y axis (m)')
                plt.show(block=False)

    def power_add(self, shape, power_type, initial_point, final_point, power):
        if shape == 'square':
            self.object.power_add(initial_point, final_point, power, shape=shape, power_type=power_type)
        if shape == 'circle':
            self.object.power_add(initial_point, final_point, power, shape=shape, power_type=power_type)
        if self.draw:
            if 'Q' in self.draw and power_type == 'Q':
                self.im_Q.set_array(self.object.Q)
                vmax = max(max(self.object.Q, key=max))
                vmin = min(min(self.object.Q, key=min))
                self.im_Q.set_clim(vmin=vmin)
                self.im_Q.set_clim(vmax=vmax)
                self.figure_Q.canvas.draw()
            if 'Q0' in self.draw and power_type == 'Q0':
                self.im_Q0.set_array(self.object.Q0)
                vmax = max(max(self.object.Q0, key=max))
                vmin = min(min(self.object.Q0, key=min))
                self.im_Q0.set_clim(vmin=vmin)
                self.im_Q0.set_clim(vmax=vmax)
                self.figure_Q0.canvas.draw()

    def compute(self, time_interval, write_interval, solver='explicit_general',
                verbose=True):
        """Compute the thermal process.

        Computes the system for timeInterval, and writes into the file_name
        file every write_interval time steps. Four different solvers can be
        used: 'explicit_general', 'explicit_k(x)', 'implicit_general',
        and 'implicit_k(x)'. heat_points is a list that defines the points
        where the heat flux are calculated if modeTemp is True the compute
        method stops when the point modeTempPoint changes numFlag relative to
        the initial value

        """
        # number of time steps for the given timeInterval
        nt = int(time_interval / self.object.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for k in range(nt):

            # updates the time_passed
            self.time_passed = self.time_passed + self.object.dt

            # defines the material properties accoring to the state list
            for i in range(self.object.size[0]):
                for j in range(self.object.size[1]):
                    if self.object.state[i][j] is True:
                        self.object.rho[i][j] = self.object.materials[self.object.materials_index[i][j]].rhoa(
                            self.object.temperature[i][j][0])
                        self.object.Cp[i][j] = self.object.materials[self.object.materials_index[i][j]].cpa(
                            self.object.temperature[i][j][0])
                        self.object.k[i][j] = self.object.materials[self.object.materials_index[i][j]].ka(
                            self.object.temperature[i][j][0])
                    if self.object.state[i][j] is False:
                        self.object.rho[i][j] = self.object.materials[self.object.materials_index[i][j]].rho0(
                            self.object.temperature[i][j][0])
                        self.object.Cp[i][j] = self.object.materials[self.object.materials_index[i][j]].cp0(
                            self.object.temperature[i][j][0])
                        self.object.k[i][j] = self.object.materials[self.object.materials_index[i][j]].k0(
                            self.object.temperature[i][j][0])

            # SOLVERS

            # implicit k constant
            if solver == 'implicit_general':
                self.object.temperature, self.object.lheat = solvers.implicit_general(self.object)

            # implicit k dependent on x
            if solver == 'implicit_k(x)':
                self.object.temperature, self.object.lheat = solvers.implicit_k(self.object)

            # explicit k constant
            if solver == 'explicit_general':
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])
                self.object.temperature, self.object.lheat = solvers.explicit_general(self.object)
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])

            # explicit k dependent on x
            if solver == 'explicit_k(x)':
                self.object.temperature, self.object.lheat = solvers.explicit_k(self.object)

            nw = nw + 1

            if self.object.file_name:
                if nw + 1 == write_interval or k == 0 or k == nt - 1:
                    line = '%f' % self.time_passed
                    for i in range(self.object.size[0]):
                        for j in range(self.object.size[1]):
                            new_line = ',%f' % self.object.temperature[i][j][1]
                            line = line + new_line
                        f = open(self.object.file_name, 'a')
                        f.write(line+'\n')
                        f.close()
            if self.draw:
                for drawing in self.draw:
                    if drawing == 'temperature':
                        if nw + 1 == write_interval or k == 0 or k == nt - 1:
                            temp = []
                            for i in range(self.object.size[0]):
                                temp.append([])
                                for j in range(self.object.size[1]):
                                    temp[-1].append(self.object.temperature[i][j][0])
                            self.im.set_array(temp)
                            if not self.draw_scale:
                                vmax = max(max(temp, key=max))
                                vmin = min(min(temp, key=min))
                                self.im.set_clim(vmin=vmin)
                                self.im.set_clim(vmax=vmax)
                            self.figure.canvas.draw()

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('pogress:', int(100*k/nt), '%', end="\r")

        if verbose:
            print('Finished simulation')
