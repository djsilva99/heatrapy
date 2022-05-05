"""Contains the class single_object.

Used to compute single two-dimensional system models.

"""

import copy
from .. import solvers
from . import Object
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


class SingleObject:
    """Single_object class.

    This class solves numerically the two-dimensional heat conduction equation.

    """

    def __init__(self, amb_temperature, material='Cu', dx=0.01, dy=0.01,
                 dt=0.1, size=(10, 10), file_name=None,
                 boundaries=(0, 0, 0, 0), Q=[], Q0=[], initial_state=False,
                 materials_path=False,
                 draw=['temperature', 'materials'], draw_scale=None):
        """Object initialization.

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
        standard heatrapy database. `draw` is a list of strings representing
        the online plots. In this version live plotting can be performed for
        `temperature`, `materials`, `state`, `Q` and `Q0`. If the list is
        empty, then no drawing is performed. `draw_scale` is a list of two
        values, representing the minimum and maximum temperature to be drawn.
        If None, there are no limits. `Q` is a list of fixed heat source
        coefficient and `Q0` is a list of temperature dependent heat source
        coefficient.

        """
        boundaries = tuple(boundaries)
        cond01 = isinstance(amb_temperature, float)
        cond01 = cond01 or isinstance(amb_temperature, int)
        cond02 = isinstance(material, str)
        cond05 = isinstance(dx, int) or isinstance(dx, float)
        cond06 = isinstance(dy, int) or isinstance(dy, float)
        cond07 = isinstance(dt, int) or isinstance(dt, float)
        cond08 = isinstance(file_name, str)
        cond08 = cond08 or (file_name is None)
        cond09 = isinstance(boundaries, tuple)
        cond10 = isinstance(initial_state, bool)
        cond11 = isinstance(Q, list)
        cond12 = isinstance(Q0, list)
        cond13 = isinstance(draw, list)
        cond14 = isinstance(draw_scale, list) or isinstance(draw_scale, tuple)
        cond14 = cond14 or draw_scale is None
        condition = cond01 and cond02 and cond05
        condition = condition and cond06 and cond07 and cond08 and cond09
        condition = condition and cond10 and cond11 and cond12 and cond13
        condition = condition and cond14
        if not condition:
            raise ValueError

        self.materials_path = materials_path
        self.time_passed = 0.
        self.size = size
        self.dt = dt
        self.dx = dx
        self.dy = dy
        self.object = Object(amb_temperature, material=material, dx=dx, dy=dy,
                             dt=dt, size=size, file_name=file_name,
                             boundaries=boundaries, Q=[], Q0=[],
                             initial_state=initial_state,
                             materials_path=materials_path)

        # Initializes plotting
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
                    data.append((1-t[0], t[2], t[1]))
                reverse.append(sorted(data))
            linear = dict(zip(k, reverse))
            my_cmap_r = matplotlib.colors.LinearSegmentedColormap(name, linear)
            cmap_r = my_cmap_r
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
                        extent = [0, size[0]*dx, 0, size[1]*dy]
                        self.im = self.ax.imshow(np.transpose(temp), vmax=vmax,
                                                 vmin=vmin, cmap=cmap_r,
                                                 extent=extent, origin='lower',
                                                 interpolation='hamming')
                    else:
                        temp = np.array(temp)
                        extent = [0, size[0]*dx, 0, size[1]*dy]
                        self.im = self.ax.imshow(np.transpose(temp),
                                                 vmax=self.draw_scale[0],
                                                 vmin=self.draw_scale[1],
                                                 cmap=cmap_r, extent=extent,
                                                 origin='lower',
                                                 interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax.figure.colorbar(self.im, ax=self.ax,
                                                   **cbar_kw)
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
                    cmap = plt.get_cmap("gray", 2)
                    extent = [0, size[0]*dx, 0, size[1]*dy]
                    self.im_state = self.ax_state.imshow(np.transpose(state),
                                                         vmax=1.5, vmin=-0.5,
                                                         cmap=cmap,
                                                         extent=extent,
                                                         origin='lower')
                    cbar_kw = {}
                    cbarlabel = ""
                    qrates = np.array(['active', 'inactive'])
                    value = np.linspace(0, 1, 2)
                    norm = matplotlib.colors.BoundaryNorm(value, 1)
                    func = lambda x, pos: qrates[::-1][norm(x)]
                    fmt = matplotlib.ticker.FuncFormatter(func)
                    cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
                    cbar_stat = self.ax_state.figure.colorbar(self.im_state,
                                                              ax=self.ax_state,
                                                              **cbar_kw)
                    cbar_stat.ax.set_ylabel(cbarlabel, rotation=-90,
                                            va="bottom")
                    self.ax_state.set_title('State')
                    self.ax_state.set_xlabel('x axis (m)')
                    self.ax_state.set_ylabel('y axis (m)')
                    plt.show(block=False)
                value = len(self.object.materials_name) > 1
                if drawing == 'materials' and value:
                    self.figure_materials = plt.figure()
                    self.ax_materials = self.figure_materials.add_subplot(111)
                    vmax = len(self.object.materials)-1
                    vmin = 0
                    material_id = np.array(self.object.materials_index)
                    cmap = plt.get_cmap("gray", vmax)
                    extent = [0, size[0]*dx, 0, size[1]*dy]
                    origin = 'lower'
                    value = np.transpose(material_id)
                    self.im_materials = self.ax_materials.imshow(value,
                                                                 vmax=vmax,
                                                                 vmin=vmin,
                                                                 cmap=cmap,
                                                                 extent=extent,
                                                                 origin=origin)
                    cbar_kw = {}
                    cbarlabel = ""
                    qrates = np.array(self.object.materials_name)
                    value_1 = np.linspace(0, 1, len(self.object.materials))
                    value_2 = len(self.object.materials)-1
                    norm = matplotlib.colors.BoundaryNorm(value_1, value_2)
                    func = lambda x, pos: qrates[::-1][norm(x)]
                    fmt = matplotlib.ticker.FuncFormatter(func)
                    cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
                    value_1 = self.im_materials
                    value_2 = self.ax_materials
                    cbar_m = self.ax_materials.figure.colorbar(value_1,
                                                               ax=value_2,
                                                               **cbar_kw)
                    cbar_m.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
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
                    extent = [0, size[0]*dx, 0, size[1]*dy]
                    self.im_Q = self.ax_Q.imshow(np.transpose(temp), vmax=vmax,
                                                 vmin=vmin, cmap='inferno',
                                                 extent=extent, origin='lower',
                                                 interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax_Q.figure.colorbar(self.im_Q, ax=self.ax_Q,
                                                     **cbar_kw)
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
                    extent = [0, size[0]*dx, 0, size[1]*dy]
                    self.im_Q0 = self.ax_Q0.imshow(np.transpose(temp),
                                                   vmax=vmax, vmin=vmin,
                                                   cmap='inferno',
                                                   extent=extent,
                                                   origin='lower',
                                                   interpolation='hamming')
                    cbar_kw = {}
                    cbar = self.ax_Q0.figure.colorbar(self.im_Q0,
                                                      ax=self.ax_Q0,
                                                      **cbar_kw)
                    self.ax_Q0.set_title('Q0 (W/m²)')
                    self.ax_Q0.set_xlabel('x axis (m)')
                    self.ax_Q0.set_ylabel('y axis (m)')
                    plt.show(block=False)

    def show_figure(self, figure_type):
        """Plotting.

        Initializes a specific plotting. `figure_type` is a string identifying
        the plotting. In this version live plotting can be performed for
        `temperature`, `materials`, `state`, `Q` and `Q0`.

        """
        # check the validity of inputs
        condition = isinstance(figure_type, str)
        if not condition:
            raise ValueError

        if figure_type not in self.draw:
            self.draw.append(figure_type)
        if figure_type == 'temperature':
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
                extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
                self.im = self.ax.imshow(np.transpose(temp), vmax=vmax,
                                         vmin=vmin, cmap='jet', extent=extent,
                                         origin='lower',
                                         interpolation='hamming')
            else:
                temp = np.array(temp)
                extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
                self.im = self.ax.imshow(np.transpose(temp),
                                         vmax=self.draw_scale[0],
                                         vmin=self.draw_scale[1],
                                         cmap='jet', extent=extent,
                                         origin='lower',
                                         interpolation='hamming')
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
            extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
            self.im_state = self.ax_state.imshow(np.transpose(state), vmax=1.5,
                                                 vmin=-0.5,
                                                 cmap=plt.get_cmap("gray", 2),
                                                 extent=extent, origin='lower')
            cbar_kw = {}
            cbarlabel = ""  # "state"
            qrates = np.array(['active', 'inactive'])
            norm = matplotlib.colors.BoundaryNorm(np.linspace(0, 1, 2), 1)
            func = lambda x, pos: qrates[::-1][norm(x)]
            fmt = matplotlib.ticker.FuncFormatter(func)
            cbar_kw = dict(ticks=np.arange(-3, 4), format=fmt)
            cbar_state = self.ax_state.figure.colorbar(self.im_state,
                                                       ax=self.ax_state,
                                                       **cbar_kw)
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
            cmap = plt.get_cmap("PiYG", vmax+1)
            extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
            value = np.transpose(material_id)
            self.im_materials = self.ax_materials.imshow(value,
                                                         vmax=vmax,
                                                         vmin=vmin,
                                                         cmap=cmap,
                                                         extent=extent,
                                                         origin='lower')
            cbar_kw = {}
            cbarlabel = ""
            materials_name_list = copy.deepcopy(self.object.materials_name)
            materials_name_list.reverse()
            qrates = np.array(materials_name_list)
            value = np.linspace(0, len(self.object.materials) - 1,
                                len(self.object.materials))
            norm = matplotlib.colors.BoundaryNorm(value,
                                                  len(self.object.materials)-1)
            func = lambda x, pos: qrates[::-1][norm(x)]
            fmt = matplotlib.ticker.FuncFormatter(func)
            cbar_kw = dict(ticks=np.arange(0, len(self.object.materials)+1),
                           format=fmt)
            cbar_m = self.ax_materials.figure.colorbar(self.im_materials,
                                                       ax=self.ax_materials,
                                                       **cbar_kw)
            cbar_m.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
            self.ax_materials.set_title('Materials')
            self.ax_materials.set_xlabel('x axis (m)')
            self.ax_materials.set_ylabel('y axis (m)')
            plt.show(block=False)
            print(self.object.materials_name)
        if figure_type == 'materials' and len(self.object.materials_name) == 1:
            print(self.object.materials_name)
        if figure_type == 'Q':
            self.figure_Q = plt.figure()
            self.ax_Q = self.figure_Q.add_subplot(111)
            temp = np.array(self.object.Q)
            vmax = 1
            vmin = 0
            extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
            self.im_Q = self.ax_Q.imshow(np.transpose(temp), vmax=vmax,
                                         vmin=vmin, cmap='inferno',
                                         extent=extent, origin='lower',
                                         interpolation='hamming')
            cbar_kw = {}
            cbar = self.ax_Q.figure.colorbar(self.im_Q, ax=self.ax_Q,
                                             **cbar_kw)
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
            extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
            self.im_Q0 = self.ax_Q0.imshow(np.transpose(temp), vmax=vmax,
                                           vmin=vmin, cmap='inferno',
                                           extent=extent,
                                           origin='lower',
                                           interpolation='hamming')
            cbar_kw = {}
            cbar = self.ax_Q0.figure.colorbar(self.im_Q0, ax=self.ax_Q0,
                                              **cbar_kw)
            self.ax_Q0.set_title('Q0 (W/m²)')
            self.ax_Q0.set_xlabel('x axis (m)')
            self.ax_Q0.set_ylabel('y axis (m)')
            plt.show(block=False)

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

        self.object.activate(initial_point=initial_point,
                             final_point=final_point, shape=shape)
        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    temp = []
                    for i in range(self.object.size[0]):
                        temp.append([])
                        for j in range(self.object.size[1]):
                            temp[-1].append(self.object.temperature[i][j][0])
                    self.im.set_array(np.transpose(temp))
                    if not self.draw_scale:
                        vmax = max(max(temp, key=max))
                        vmin = min(min(temp, key=min))
                        self.im.set_clim(vmin=vmin)
                        self.im.set_clim(vmax=vmax)
                    self.figure.canvas.draw()
                if drawing == 'state':
                    self.im_state.set_array(np.transpose(self.object.state))
                    self.figure_state.canvas.draw()

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

        self.object.deactivate(initial_point=initial_point,
                               final_point=final_point, shape=shape)
        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    temp = []
                    for i in range(self.object.size[0]):
                        temp.append([])
                        for j in range(self.object.size[1]):
                            temp[-1].append(self.object.temperature[i][j][0])
                    self.im.set_array(np.transpose(temp))
                    if not self.draw_scale:
                        vmax = max(max(temp, key=max))
                        vmin = min(min(temp, key=min))
                        self.im.set_clim(vmin=vmin)
                        self.im.set_clim(vmax=vmax)
                    self.figure.canvas.draw()
                if drawing == 'state':
                    self.im_state.set_array(np.transpose(self.object.state))
                    self.figure_state.canvas.draw()

    def change_boundaries(self, boundaries):
        """Boundary change.

        Changes the `boundaries` variable.

        """
        # check the validity of inputs
        if isinstance(boundaries, tuple):
            if len(boundaries) == 4:
                condition = True
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        self.object.boundaries = boundaries

    def change_material(self, shape, material, initial_point, length,
                        state=False):
        """Material change.

        Changes the material of a given piece of the background material. If
        `shape` is `'square'`, then the initial_point is the tuple (x,y) of the
        bottom left point and the `length` is the tuple (x,y) of the length. If
        the shape is `'circle'`, the `initial_point` is the tuple (x,y) of the
        center of the circle and `length` is its radius.

        """
        # check the validity of inputs
        if isinstance(shape, str):
            if shape == 'square':
                value = isinstance(initial_point, tuple)
                if value and isinstance(length, tuple):
                    condition = len(initial_point) == 2
                    condition = condition and len(length) == 2
                else:
                    condition = False
            elif shape == 'circle':
                value = isinstance(length, int)
                value = value or isinstance(length, float)
                value = value and isinstance(initial_point, tuple)
                if value:
                    condition = len(initial_point) == 2
                else:
                    condition = False
            else:
                condition = False
        else:
            condition = False
        condition = condition and isinstance(state, bool)
        if not condition:
            raise ValueError

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
        except:
            pass

        if self.draw:
            if 'materials' in self.draw:
                self.figure_materials = plt.figure()
                self.ax_materials = self.figure_materials.add_subplot(111)
                vmax = len(self.object.materials)-1
                vmin = 0
                cmap = plt.get_cmap("PiYG", vmax+1)
                extent = [0, self.size[0]*self.dx, 0, self.size[1]*self.dy]
                material_id = np.array(self.object.materials_index)
                value = np.transpose(material_id)
                self.im_materials = self.ax_materials.imshow(value,
                                                             vmax=vmax,
                                                             vmin=vmin,
                                                             cmap=cmap,
                                                             extent=extent,
                                                             origin='lower')
                cbar_kw = {}
                cbarlabel = ""
                materials_name_list = copy.deepcopy(self.object.materials_name)
                materials_name_list.reverse()
                qrates = np.array(materials_name_list)
                value_1 = np.linspace(0, len(self.object.materials)-1,
                                      len(self.object.materials))
                value_2 = len(self.object.materials)-1
                norm = matplotlib.colors.BoundaryNorm(value_1, value_2)
                func = lambda x, pos: qrates[::-1][norm(x)]
                fmt = matplotlib.ticker.FuncFormatter(func)
                value = np.arange(0, len(self.object.materials)+1)
                cbar_kw = dict(ticks=value, format=fmt)
                cba_m = self.ax_materials.figure.colorbar(self.im_materials,
                                                          ax=self.ax_materials,
                                                          **cbar_kw)
                cba_m.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                self.ax_materials.set_title('Materials')
                self.ax_materials.set_xlabel('x axis (m)')
                self.ax_materials.set_ylabel('y axis (m)')
                plt.show(block=False)

    def change_power(self, shape, power_type, initial_point, final_point,
                     power):
        """Power change.

        Changes the power matrix of the thermal object. If `shape` is
        `'square'`, then the `initial_point` is the tuple (x,y) of the bottom
        left point and the `final_point` is the tuple (x,y) of the top right
        point. If the `shape` is `'circle'`, the `initial_point` is the tuple
        (x,y) of the center of the circle and `final_point` is its radius.
        `power` is the value of the power to add, and `power_type` is the type
        of power to be introduced, which has the value `'Q'` if it is
        temperature dependent and `'Q0'` if it is temperature independent.

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

        if shape == 'square':
            self.object.power_add(initial_point, final_point, power,
                                  shape=shape, power_type=power_type)
        if shape == 'circle':
            self.object.power_add(initial_point, final_point, power,
                                  shape=shape, power_type=power_type)
        if self.draw:
            if 'Q' in self.draw and power_type == 'Q':
                self.im_Q.set_array(np.transpose(self.object.Q))
                vmax = max(max(self.object.Q, key=max))
                vmin = min(min(self.object.Q, key=min))
                self.im_Q.set_clim(vmin=vmin)
                self.im_Q.set_clim(vmax=vmax)
                self.figure_Q.canvas.draw()
            if 'Q0' in self.draw and power_type == 'Q0':
                self.im_Q0.set_array(np.transpose(self.object.Q0))
                vmax = max(max(self.object.Q0, key=max))
                vmin = min(min(self.object.Q0, key=min))
                self.im_Q0.set_clim(vmin=vmin)
                self.im_Q0.set_clim(vmax=vmax)
                self.figure_Q0.canvas.draw()

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
        all_solvers = ['explicit_general', 'explicit_k(x)']
        cond3 = isinstance(solver, str) and solver in all_solvers
        cond4 = isinstance(verbose, bool)
        condition = cond1 and cond2 and cond3 and cond4
        if not condition:
            raise ValueError

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
                        ix = self.object.materials_index[i][j]
                        self.object.rho[i][j] = self.object.materials[ix].rhoa(
                            self.object.temperature[i][j][0])
                        self.object.Cp[i][j] = self.object.materials[ix].cpa(
                            self.object.temperature[i][j][0])
                        self.object.k[i][j] = self.object.materials[ix].ka(
                            self.object.temperature[i][j][0])
                    if self.object.state[i][j] is False:
                        ix = self.object.materials_index[i][j]
                        self.object.rho[i][j] = self.object.materials[ix].rho0(
                            self.object.temperature[i][j][0])
                        self.object.Cp[i][j] = self.object.materials[ix].cp0(
                            self.object.temperature[i][j][0])
                        self.object.k[i][j] = self.object.materials[ix].k0(
                            self.object.temperature[i][j][0])

            # SOLVERS
            # explicit k constant
            if solver == 'explicit_general':
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])
                value = solvers.explicit_general(self.object)
                self.object.temperature, self.object.lheat = value
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])
            # explicit k constant
            if solver == 'explicit_k(x)':
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])
                value = solvers.explicit_k(self.object)
                self.object.temperature, self.object.lheat = value
                temp = []
                for i in range(self.object.size[0]):
                    temp.append([])
                    for j in range(self.object.size[1]):
                        temp[-1].append(self.object.temperature[i][j][0])

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
                                    value = self.object.temperature[i][j][0]
                                    temp[-1].append(value)
                            try:
                                self.im.set_array(np.transpose(temp))
                                if not self.draw_scale:
                                    vmax = max(max(temp, key=max))
                                    vmin = min(min(temp, key=min))
                                    self.im.set_clim(vmin=vmin)
                                    self.im.set_clim(vmax=vmax)
                                self.figure.canvas.draw()
                            except:
                                pass

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('pogress:', int(100*k/nt), '%', end="\r")
        if verbose:
            print('Finished simulation')
