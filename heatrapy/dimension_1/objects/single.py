"""Contains the class single_object.

Used to compute single thermal objects.

"""

from .. import solvers
from . import Object
import matplotlib.pyplot as plt
import numpy as np


class SingleObject:
    """Single_object class.

    This class solves numerically the heat conduction equation for 1 dimension
    of a single material(s). The class has 6 methods.

    """

    def __init__(self, amb_temperature, materials=('Cu',), borders=(1, 11),
                 materials_order=(0,), dx=0.01, dt=0.1, file_name=None,
                 boundaries=(0, 0), initial_state=False,
                 materials_path=False, draw=['temperature'], draw_scale=None):
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
        heatrapy database. `draw` is a list of strings representing the online
        plots. In this version only `'temperature'` can be potted. If the list
        is empty, then no drawing is performed. `draw_scale` is a list of two
        values, representing the minimum and maximum temperature to be drawn.
        If None, there are no limits.

        """
        # check the validity of inputs
        materials = tuple(materials)
        borders = tuple(borders)
        materials_order = tuple(materials_order)
        boundaries = tuple(boundaries)
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
        cond10 = isinstance(initial_state, bool)
        if isinstance(draw, list):
            cond15 = True
        elif draw is None:
            cond15 = True
        else:
            cond15 = False
        if isinstance(draw_scale, list) or isinstance(draw_scale, tuple):
            cond16 = (len(draw_scale) == 2)
        elif draw_scale is None:
            cond16 = True
        else:
            cond16 = False
        condition = cond01 and cond02 and cond03 and cond04 and cond05
        condition = condition and cond06 and cond07 and cond08
        condition = condition and cond10
        condition = condition and cond15 and cond16
        if not condition:
            raise ValueError

        self.object = Object(amb_temperature, materials=materials,
                             borders=borders, materials_order=materials_order,
                             dx=dx, dt=dt, file_name=file_name,
                             boundaries=boundaries,
                             initial_state=initial_state,
                             materials_path=materials_path)

        # initializes the plotting
        self.draw = draw
        self.draw_scale = draw_scale
        for drawing in self.draw:
            if drawing == 'temperature':
                self.figure = plt.figure()
                self.ax = self.figure.add_subplot(111)
                temp = []
                for i in range(len(self.object.temperature)):
                    temp.append(self.object.temperature[i][0])
                if not self.draw_scale:
                    vmax = max(temp)
                    vmin = min(temp)
                    if vmax == vmin:
                        vmin = vmin - 0.1
                        vmax = vmax + 0.1
                    temp = np.array(temp)
                    x_plot = [self.object.dx*j for j in range(len(temp))]
                    self.online, = self.ax.plot(x_plot, temp)
                    self.ax.set_ylim([vmin, vmax])
                else:
                    temp = np.array(temp)
                    x_plot = [self.object.dx*j for j in range(len(temp))]
                    self.online, = self.ax.plot(x_plot, temp)
                    self.ax.set_ylim(self.draw_scale)
                self.ax.set_title('Temperature (K)')
                self.ax.set_xlabel('x axis (m)')
                self.ax.set_ylabel('temperature (K)')
                plt.show(block=False)

    def show_figure(self, figure_type, draw_scale=None):
        """Plotting.

        Initializes a specific live plotting. `figure_type` is a string
        identifying the plotting. This version only allows the plotting of the
        'temperature'. `draw_scale` defines the range of temperatures. If None,
        this range is found automatically for every frame.

        """
        # check the validity of inputs
        if isinstance(draw_scale, list) or isinstance(draw_scale, tuple):
            condition = (len(draw_scale) == 2)
        elif draw_scale is None:
            condition = True
        else:
            condition = False
        condition = condition and isinstance(figure_type, str)
        if not condition:
            raise ValueError

        self.draw_scale = draw_scale
        if figure_type == 'temperature':
            if figure_type not in self.draw:
                self.draw.append(figure_type)
            self.figure = plt.figure()
            self.ax = self.figure.add_subplot(111)
            temp = []
            for i in range(len(self.object.temperature)):
                temp.append(self.object.temperature[i][0])
            if not self.draw_scale:
                vmax = max(temp)
                vmin = min(temp)
                if vmax == vmin:
                    vmin = vmin - 0.1
                    vmax = vmax + 0.1
                temp = np.array(temp)
                x_plot = [self.object.dx*j for j in range(len(temp))]
                self.online, = self.ax.plot(x_plot, temp)
                self.ax.set_ylim([vmin, vmax])
            else:
                temp = np.array(temp)
                x_plot = [self.object.dx*j for j in range(len(temp))]
                self.online, = self.ax.plot(x_plot, temp)
                self.ax.set_ylim(self.draw_scale)
            self.ax.set_title('Temperature (K)')
            self.ax.set_xlabel('x axis (m)')
            self.ax.set_ylabel('temperature (K)')
            plt.show(block=False)

    def activate(self, initial_point, final_point):
        """Activation.

        Activates the thermal object between `initial_point` to `final_point`.
        """
        # check the validity of inputs
        condition = isinstance(initial_point, int)
        condition = condition and isinstance(final_point, int)
        if not condition:
            raise ValueError

        self.object.activate(initial_point, final_point)

        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    try:
                        temp = []
                        for i in range(len(self.object.temperature)):
                            temp.append(self.object.temperature[i][0])
                        if not self.draw_scale:
                            vmax = max(temp)
                            vmin = min(temp)
                            if vmax == vmin:
                                vmin = vmin - 0.1
                                vmax = vmax + 0.1
                            temp = np.array(temp)
                            self.online.set_ydata(temp)
                            self.ax.set_ylim([vmin, vmax])
                        else:
                            temp = np.array(temp)
                            self.online.set_ydata(temp)
                        self.figure.canvas.draw()
                    except:
                        pass

    def deactivate(self, initial_point, final_point):
        """Deactivation.

        Deactivates the thermal object between `initial_point` to
        `final_point`.
        """
        # check the validity of inputs
        condition = isinstance(initial_point, int)
        condition = condition and isinstance(final_point, int)
        if not condition:
            raise ValueError

        self.object.deactivate(initial_point, final_point)

        if self.draw:
            for drawing in self.draw:
                if drawing == 'temperature':
                    try:
                        temp = []
                        for i in range(len(self.object.temperature)):
                            temp.append(self.object.temperature[i][0])
                        if not self.draw_scale:
                            vmax = max(temp)
                            vmin = min(temp)
                            if vmax == vmin:
                                vmin = vmin - 0.1
                                vmax = vmax + 0.1
                            temp = np.array(temp)
                            self.online.set_ydata(temp)
                            self.ax.set_ylim([vmin, vmax])
                        else:
                            temp = np.array(temp)
                            self.online.set_ydata(temp)
                        self.figure.canvas.draw()
                    except:
                        pass

    def change_power(self, power_type, power, initial_point, final_point):
        """Heat power source change.

        Changes the coeficients for the heat power sources by a value of power
        from `initial_point` to `final_point`. `power_type` is a string that
        represents the type of coefficient, i.e. 'Q' or 'Q0'.

        """
        # check the validity of inputs
        value = isinstance(initial_point, int)
        if value and isinstance(final_point, int):
            cond1 = True
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
        if not (cond1 and cond2 and cond3):
            raise ValueError

        if power_type == 'Q':
            for j in range(initial_point, final_point):
                self.object.Q[j] = power
        if power_type == 'Q0':
            for j in range(initial_point, final_point):
                self.object.Q0[j] = power

    def change_boundaries(self, boundaries):
        """Boundary change.

        Changes the `boundaries` variable.

        """
        # check the validity of inputs
        if isinstance(boundaries, tuple):
            if len(boundaries) == 2:
                condition = True
            else:
                condition = False
        else:
            condition = False
        if not condition:
            raise ValueError

        self.object.boundaries = boundaries

    def compute(self, time_interval, write_interval, solver='explicit_k(x)',
                verbose=True):
        """Compute the thermal process.

        Computes the system for time_interval seconds, and writes into the
        `file_name` file every `write_interval` time steps. Four different
        solvers can be used: `'explicit_general'`, `'explicit_k(x)'`,
        `'implicit_general'`, and `'implicit_k(x)'`. If `verbose = True`, then
        the progress of the computation progress is shown.

        """
        # check the validity of inputs
        cond1 = isinstance(time_interval, float)
        cond1 = cond1 or isinstance(time_interval, int)
        cond2 = isinstance(write_interval, int)
        if isinstance(solver, str):
            all_solvers = ['implicit_general', 'implicit_k(x)',
                           'explicit_k(x)', 'explicit_general']
            if solver in all_solvers:
                cond3 = True
            else:
                cond3 = False
        else:
            cond3 = False
        cond4 = isinstance(verbose, bool)
        condition = cond1 and cond2 and cond3 and cond4
        if not condition:
            raise ValueError

        # number of time steps for the given timeInterval
        nt = int(time_interval / self.object.dt)

        # number of time steps counting from the last writing process
        nw = 0

        # computes
        for j in range(nt):

            # updates the time_passed
            self.object.time_passed = self.object.time_passed + self.object.dt

            # defines the material properties accoring to the state list
            for i in range(1, self.object.num_points - 1):
                if self.object.state[i] is True:
                    value = self.object.materials_index[i]
                    self.object.rho[i] = self.object.materials[value].rhoa(
                        self.object.temperature[i][0])
                    self.object.Cp[i] = self.object.materials[value].cpa(
                        self.object.temperature[i][0])
                    self.object.k[i] = self.object.materials[value].ka(
                        self.object.temperature[i][0])
                if self.object.state[i] is False:
                    value = self.object.materials_index[i]
                    self.object.rho[i] = self.object.materials[value].rho0(
                        self.object.temperature[i][0])
                    self.object.Cp[i] = self.object.materials[value].cp0(
                        self.object.temperature[i][0])
                    self.object.k[i] = self.object.materials[value].k0(
                        self.object.temperature[i][0])

            # SOLVERS
            # implicit k constant
            if solver == 'implicit_general':
                value = solvers.implicit_general(self.object)
                self.object.temperature, self.object.lheat = value
            # implicit k dependent on x
            if solver == 'implicit_k(x)':
                value = solvers.implicit_k(self.object)
                self.object.temperature, self.object.lheat = value
            # explicit k constant
            if solver == 'explicit_general':
                value = solvers.explicit_general(self.object)
                self.object.temperature, self.object.lheat = value
            # explicit k dependent on x
            if solver == 'explicit_k(x)':
                value = solvers.explicit_k(self.object)
                self.object.temperature, self.object.lheat = value

            nw = nw + 1

            if self.draw:
                for drawing in self.draw:
                    if drawing == 'temperature':
                        try:
                            value = nw + 1 == write_interval
                            if value or j == 0 or j == nt - 1:
                                temp = []
                                for i in range(len(self.object.temperature)):
                                    temp.append(self.object.temperature[i][0])
                                if not self.draw_scale:
                                    vmax = max(temp)
                                    vmin = min(temp)
                                    if vmax == vmin:
                                        vmin = vmin - 0.1
                                        vmax = vmax + 0.1
                                    temp = np.array(temp)
                                    self.online.set_ydata(temp)
                                    self.ax.set_ylim([vmin, vmax])
                                else:
                                    temp = np.array(temp)
                                    self.online.set_ydata(temp)
                                self.figure.canvas.draw()
                        except:
                            pass

            # writes the temperature to file_name file ...
            # if the number of time steps is verified
            if self.object.file_name:
                if nw == write_interval or j == 0 or j == nt - 1:
                    line = '%f,' % self.object.time_passed
                    for i in self.object.temperature:
                        new_line = '%f,' % i[1]
                        line = line + new_line
                    line = line[:-1] + '\n'
                    f = open(self.object.file_name, 'a')
                    f.write(line)
                    f.close()

            if nw == write_interval:
                nw = 0
                if verbose:
                    print('pogress:', int(100*j/nt), '%', end="\r")

        if verbose:
            print('Finished simulation')
