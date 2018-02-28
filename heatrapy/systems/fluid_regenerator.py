# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class magcalsys_fluidAMR_1D.

Used to compute 1-dimensional models for ferroic-based systems using fluids
with the regenerative heat process.

"""

from .. import objects
import time
import numpy as np
from .. import fields


class magcalsys_fluidAMR_1D:

    """magcalsys_fluidAMR_1D class

    Computes the thermal processes for 1-dimensional ferroic-based systems
    using fluids. The active regenerative processes can be used with the
    several allowed modes for application and removal of fields. Cascades of
    materials can also be computed.

    """

    def __init__(self, fileName, ambTemperature=293, fluid_length=100,
                 MCM_length=20, rightReservoir_length=3,
                 leftReservoir_length=3, MCM_material=((0.000, 'Gd'),),
                 fluid_material='water', leftReservoir_material='Cu',
                 rightReservoir_material='Cu', freq=1, dt=0.001, dx=0.004,
                 stopCriteria=1e-3, solverMode='implicit_k(x)',
                 minCycleNumber=1, maxCycleNumber=50, cyclePoints=25,
                 note=None, boundaries=((2, 293), (3, 293)), mode='heat_pump',
                 version=None, leftHEXpositions=15, rightHEXpositions=15,
                 startingField='magnetization', temperatureSensor=[3, -3],
                 demagnetizationSteps=1, magnetizationSteps=1,
                 demagnetizationMode='accelerated_left',
                 magnetizationMode='accelerated_right',
                 applied_static_field_time_ratio=(0., 0.),
                 removed_static_field_time_ratio=(0., 0.),
                 h_mcm_fluid=5000000, h_leftreservoir_fluid=5000000,
                 h_rightreservoir_fluid=5000000, mcm_discontinuity='default',
                 type_study='fixed temperature span', velocity=.2,
                 mod_freq='default'):
        """Initialization.

        fileName: file name where the temperature and heat flux are saved
        ambTemperature: ambient temperature of the whole system
        fluid_length: length of the fluid
        fluid_material: string for the material of the fluid
        MCM_length: length of the magnetocaloric material
        MCM_material: string for the material of the magnetocaloric material
        leftReservoir_length: length of the left reservoir
        rightReservoir_length: length of the right reservoir
        leftReservoir_material: string for the material of the left reservoir
        rightReservoir_material: string for the material of the right reservoir
        freq: operating frequency
        dt: times step
        dx: space step
        stopCriteria: error threshold to stop the simulation
        solverMode: the solver
        minCycleNumber: minimum number of cycles that has to be computed
        maxCycleNumber: maximum number of cycles that has to be computed
        demagnetizationSteps: number of steps during the demagnetization
        magnetizationSteps: number of steps during the magnetization
        demagnetizationMode: mode of demagnetization
            modes can be constant_right, constant_left, accelerated_right,
            accelerated_left, decelerated_right, and decelerated_left
        magnetizationMode is the mode of magnetization
            modes can be constant_right, constant_left, accelerated_right,
            accelerated_left, decelerated_right, and decelerated_left
        cyclePoints: number of points recorded for each position for each cycle
        boundaries: tuple of two entries that define the boundary condition
            for tempreture. The first corresponds to the thermal obect while
            the second defines the temperature. If 0 the boundary condition is
            insulation
        temperatureSensor: list of two space indexes used to determine the
            temperature span at the end of the simulation. The first term is
            the sensor at the hot end and the second at the cold end
        mode: mode used for the power calculations (e.g. COP) performed at
            the end of the simulation
        version: heatrapy version (default is None)
        type_study: 'no load' or 'fixed temperature span'
        h_mcm_fluid: heat transfer coefficient for fluid - MCM
        h_leftreservoir_fluid: heat transfer coefficient for
            fluid - left reservoir
        h_rigthreservoir_fluid: heat transfer coefficient for
            fluid - right reservoir
        mod_freq: if not 'default', allows to modulate the frequency according
            to a specific temperature. tuple with first element as filename,
            and second the sensor point.
        velocity: velocity of the fluid
        leftHEXpositions: distance in points from the left reservoir to the MCM
        rightHEXpositions: distance in points from the right reservoir to the
            MCM.
        applied_static_field_time_ratio: tuple with the resting time ratios
            before and after the application of the field.
        removed_static_field_time_ratio: tuple with the resting time ratios
            before and after the removal of the field.
        mcm_discontinuity: if not 'default' the MCM is devided into n pieces.
            the input is a tuple where the first entry is the number of
            discontinuities, while the second is the thickness of each
            discontinuity in meters.

        """

        if startingField == 'demagnetization':
            initialState = True
        else:
            initialState = False

        cycle_number = 0
        AMR = objects.system_objects(number_objects=4,
                                     materials=[fluid_material,
                                                MCM_material[0][1],
                                                leftReservoir_material,
                                                rightReservoir_material],
                                     objects_length=[fluid_length, MCM_length,
                                                     leftReservoir_length,
                                                     rightReservoir_length],
                                     ambTemperature=293, dx=dx, dt=dt,
                                     fileName=fileName,
                                     initialState=initialState,
                                     boundaries=boundaries)

        k = 1
        for i in range(len(MCM_material)):
            from .. import mats
            import os
            tadi = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'tadi.txt'
            tadd = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'tadd.txt'
            cpa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'cpa.txt'
            cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'cp0.txt'
            k0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'k0.txt'
            ka = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'ka.txt'
            rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'rho0.txt'
            rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/' + MCM_material[i][1] + '/' + 'rhoa.txt'
            AMR.objects[1].materials.append(mats.calmatpro(tadi, tadd, cpa,
                                                           cp0, k0, ka, rho0,
                                                           rhoa))
            for j in range(k, k+int(MCM_material[i][0]/dx)):
                len_mat = len(AMR.objects[1].materials)
                AMR.objects[1].materialsIndex[j] = len_mat - 1
            k = k + int(MCM_material[i][0]/dx)

        if mcm_discontinuity != 'default':
            from .. import mats
            import os
            tadi = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/tadi.txt'
            tadd = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/tadd.txt'
            cpa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/cpa.txt'
            cp0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/cp0.txt'
            k0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/k0.txt'
            ka = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/ka.txt'
            rho0 = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/rho0.txt'
            rhoa = os.path.dirname(os.path.realpath(__file__)) + \
                '/../database/vaccum/rhoa.txt'
            AMR.objects[1].materials.append(mats.calmatpro(tadi, tadd, cpa,
                                                           cp0, k0, ka, rho0,
                                                           rhoa))

            j = 1
            for i in range(1, len(AMR.objects[1].temperature)-1):
                i_cond1 = (j * len(AMR.objects[1].temperature) /
                           (mcm_discontinuity[0] + 1) +
                           int(mcm_discontinuity[1] / (2 * dx)))
                i_cond2 = (j * len(AMR.objects[1].temperature) /
                           (mcm_discontinuity[0] + 1) -
                           int(mcm_discontinuity[1] / (2 * dx)))
                if i < j_cond1 and i >= i_cond2:
                    val = len(AMR.objects[1].materials)-1
                    AMR.objects[1].materialsIndex[i] = val
                    i_cond3 = (j * len(AMR.objects[1].temperature) /
                               (mcm_discontinuity[0] + 1) +
                               int(mcm_discontinuity[1] / (2 * dx)))
                if i == i_cond3 and j < mcm_discontinuity[0]:
                    j = j + 1

        write_interval = cyclePoints/2
        steps = int(velocity / (2 * freq * dx))
        time_step = (1 / (2. * freq)) / steps

        if int(time_step / dt) == 0:
            print 'dt or frequency too low'

        if write_interval > int(time_step / dt):
            write_interval = 1

        # information for the log file
        print ''
        print ''
        print '######################################################'
        print ''
        print '------------------------------------------------------'
        print fileName
        print '------------------------------------------------------'
        print ''
        print 'heatconpy version:', version
        print 'Module: magcalsys_fluidAMR_1D'
        if note is not None:
            print ''
            print 'Note:', note
        print ''
        print 'Mode:', mode
        print 'Fluid:', fluid_material + ' (' + str(fluid_length*dx) + ')'
        print 'MCM material:', MCM_material
        string = ' (' + str(leftReservoir_length * dx) + ')'
        print 'Left reservoir:', leftReservoir_material + string
        string = ' (' + str(rightReservoir_length * dx) + ')'
        print 'Right reservoir:', rightReservoir_material + string
        print 'Distance between MCM and left HEX:', leftHEXpositions*dx
        print 'Distance between MCM and right HEX:', rightHEXpositions*dx
        print 'dx (m):', dx
        print 'dt (s):', dt
        print 'Frequency (Hz):', freq
        if mcm_discontinuity == 'default':
            print 'Discontinuity: None'
        else:
            print 'Discontinuity:', mcm_discontinuity
        print 'MCM - fluid heat transfer coefficient:', h_mcm_fluid
        val = h_leftreservoir_fluid
        print 'Left reservoir - fluid heat transfer coefficient:', val
        val = h_rightreservoir_fluid
        print 'Right reservoir - fluid heat transfer coefficient:', val
        val = applied_static_field_time_ratio
        print 'Applied field static time ratios:', val
        val = removed_static_field_time_ratio
        print 'Removed field static time ratios:', val
        print 'Solver:', solverMode
        print 'Magnetization mode:', magnetizationMode
        print 'Magnetization steps:', magnetizationSteps
        print 'Demagnetization mode:', demagnetizationMode
        print 'Demagnetization steps:', demagnetizationSteps
        print 'Starting Field:', startingField
        print 'Boundaries:', boundaries
        print 'Ambient temperature (K):', ambTemperature
        print 'Stop criteria:', stopCriteria
        print 'Time:', time.ctime()
        print ''

        startTime = time.time()

        if type_study == 'no load':

            value1 = AMR.objects[0].temperature[temperatureSensor[0]][0]
            value2 = AMR.objects[0].temperature[temperatureSensor[1]][0]

            if startingField == 'magnetization':

                condition = True
                while condition:

                    val = AMR.objects[0].temperature[temperatureSensor[0]][0]
                    value1 = val
                    val = AMR.objects[0].temperature[temperatureSensor[1]][0]
                    value2 = val

                    if mod_freq != 'default':
                        val = AMR.objects[1].temperature[mod_freq[1]][0]
                        temperature_sensor = val
                        input = open(mod_freq[0], 'r')
                        s = input.readlines()
                        xMod = []
                        yMod = []
                        for line in s:
                            pair = line.split(',')
                            xMod.append(float(pair[0]))
                            yMod.append(float(pair[1]))
                        input.close()
                        freq = np.interp(temperature_sensor, xMod, yMod)
                        print freq, AMR.objects[1].temperature[mod_freq[1]][0]
                        steps = int(velocity / (2 * freq * dx))
                        time_step = (1 / (2. * freq)) / steps
                        if int(time_step / dt) == 0:
                            print 'dt or frequency too low'
                        if write_interval > int(time_step / dt):
                            write_interval = 1

                    # MAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if magnetizationMode in com:
                        j = magnetizationSteps - 1
                        j_cond = magnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 - (leftReservoir_length +
                                    leftHEXpositions + MCM_length +
                                    rightHEXpositions +
                                    rightReservoir_length) / 2 - steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+n+init_pos), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos+n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos+n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p+1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length/2 - (leftReservoir_length +
                                                      leftHEXpositions +
                                                      MCM_length +
                                                      rightHEXpositions +
                                                      rightReservoir_length) /
                                    2 - steps/2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos+n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1/(2.*freq) -
                                 applied_static_field_time_ratio[1]/2)
                        if n*time_step < applied_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                md = magnetizationMode
                                ap1 = applied_static_field_time_ratio[0]
                                ap2 = applied_static_field_time_ratio[1]
                                ms = magnetizationSteps
                                time_interval = fields.operating_mode(md, ap1,
                                                                      ap2, ms,
                                                                      freq, j)
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val and flag:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)

                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if magnetizationMode not in com:
                                        if j < magnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    # DEMAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if demagnetizationMode in com:
                        j = demagnetizationSteps - 1
                        j_cond = demagnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+init_pos-n), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos-n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] /
                                             (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos-n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p + 1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos-n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1 / (2. * freq) -
                                 removed_static_field_time_ratio[1] / 2)
                        if n*time_step < removed_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                v1 = demagnetizationMode
                                v2 = removed_static_field_time_ratio[0]
                                v3 = removed_static_field_time_ratio[1]
                                v4 = demagnetizationSteps
                                val = fields.operating_mode(v1, v2, v3,
                                                            v4, freq,
                                                            j)
                                time_interval = val
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = demagnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step - previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    ds = demagnetizationSteps
                                    if time_passed == 0. and j < ds and flag:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if demagnetizationMode not in com:
                                        if j < demagnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    if value1 == value2:
                        condition = cycle_number < maxCycleNumber - 1
                    else:
                        cond1 = cycle_number < minCycleNumber
                        vl1 = leftReservoir_length/2
                        vl2 = rightReservoir_length/2
                        cond2 = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                          AMR.objects[3].temperature[vl2][0]) -
                                      abs(value1-value2)))/abs(value1-value2))
                        cond2 = cond2 > stopCriteria
                        cond3 = cycle_number < maxCycleNumber - 1
                        condition = (cond1 or cond2) and cond3

                    cycle_number = cycle_number + 1

            elif startingField == 'demagnetization':
                condition = True
                while condition:
                    val = AMR.objects[0].temperature[temperatureSensor[0]][0]
                    value1 = val
                    val = AMR.objects[0].temperature[temperatureSensor[1]][0]
                    value2 = val

                    if mod_freq != 'default':
                        val = AMR.objects[1].temperature[mod_freq[1]][0]
                        temperature_sensor = val
                        input = open(mod_freq[0], 'r')
                        s = input.readlines()
                        xMod = []
                        yMod = []
                        for line in s:
                            pair = line.split(',')
                            xMod.append(float(pair[0]))
                            yMod.append(float(pair[1]))
                        input.close()
                        freq = np.interp(temperature_sensor, xMod, yMod)
                        steps = int(velocity / (2 * freq * dx))
                        time_step = (1 / (2. * freq)) / steps
                        if int(time_step / dt) == 0:
                            print 'dt or frequency too low'
                        if write_interval > int(time_step / dt):
                            write_interval = 1

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if demagnetizationMode in com:
                        j = demagnetizationSteps - 1
                        j_cond = demagnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+init_pos-n), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) /
                                        2 + steps / 2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos-n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 +
                                        steps / 2)
                            for i in range(MCM_length):
                                val1 = (p * len(AMR.objects[1].temperature) /
                                        (mcm_discontinuity[0] + 1) +
                                        int(mcm_discontinuity[1] / (2 * dx)))
                                val2 = (p * len(AMR.objects[1].temperature) /
                                        (mcm_discontinuity[0] + 1) -
                                        int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < val1 and i+1 >= val2):
                                    contacts.add(((0, i+init_pos-n), (1, i+1),
                                                  h_mcm_fluid))
                                val = (p * len(AMR.objects[1].temperature) /
                                       (mcm_discontinuity[0] + 1) +
                                       int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == val and p < mcm_discontinuity[0]:
                                    p = p+1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos-n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts
                        val = (1 / (2. * freq) -
                               removed_static_field_time_ratio[1] / 2)
                        if n*time_step < removed_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > val:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                v1 = demagnetizationMode
                                v2 = removed_static_field_time_ratio[0]
                                v3 = removed_static_field_time_ratio[1]
                                v4 = demagnetizationSteps
                                val = fields.operating_mode(v1, v2, v3,
                                                            v4, freq,
                                                            j)
                                time_interval = val
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = demagnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j+j_inc
                                        delta_t = time_step - previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    val = demagnetizationSteps
                                    if time_passed == 0. and j < val and flag:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if demagnetizationMode not in com:
                                        if j < demagnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    # MAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if magnetizationMode in com:
                        j = magnetizationSteps - 1
                        j_cond = magnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 - (leftReservoir_length +
                                    leftHEXpositions + MCM_length +
                                    rightHEXpositions +
                                    rightReservoir_length) / 2 - steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+n+init_pos), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos+n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos+n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p+1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length/2 - (leftReservoir_length +
                                                      leftHEXpositions +
                                                      MCM_length +
                                                      rightHEXpositions +
                                                      rightReservoir_length) /
                                    2 - steps/2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos+n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1/(2.*freq) -
                                 applied_static_field_time_ratio[1]/2)
                        if n*time_step < applied_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                md = magnetizationMode
                                ap1 = applied_static_field_time_ratio[0]
                                ap2 = applied_static_field_time_ratio[1]
                                ms = magnetizationSteps
                                time_interval = fields.operating_mode(md, ap1,
                                                                      ap2, ms,
                                                                      freq, j)
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val and flag:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)

                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if magnetizationMode not in com:
                                        if j < magnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    if value1 == value2:
                        condition = cycle_number < maxCycleNumber - 1
                    else:
                        cond1 = cycle_number < minCycleNumber
                        vl1 = leftReservoir_length/2
                        vl2 = rightReservoir_length/2
                        cond2 = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                          AMR.objects[3].temperature[vl2][0]) -
                                      abs(value1-value2)))/abs(value1-value2))
                        cond2 = cond2 > stopCriteria
                        cond3 = cycle_number < maxCycleNumber - 1
                        condition = (cond1 or cond2) and cond3

                    cycle_number = cycle_number + 1

            else:
                print 'incorrect startingField'

        if type_study == 'fixed temperature span':

            if startingField == 'magnetization':

                condition = True
                q1 = AMR.q1

                while condition:

                    value1 = AMR.q1 - q1
                    q1 = AMR.q1
                    q2 = AMR.q2

                    if mod_freq != 'default':
                        val = AMR.objects[1].temperature[mod_freq[1]][0]
                        temperature_sensor = val
                        input = open(mod_freq[0], 'r')
                        s = input.readlines()
                        xMod = []
                        yMod = []
                        for line in s:
                            pair = line.split(',')
                            xMod.append(float(pair[0]))
                            yMod.append(float(pair[1]))
                        input.close()
                        freq = np.interp(temperature_sensor, xMod, yMod)
                        steps = int(velocity / (2 * freq * dx))
                        time_step = (1 / (2. * freq)) / steps
                        if int(time_step / dt) == 0:
                            print 'dt or frequency too low'
                        if write_interval > int(time_step / dt):
                            write_interval = 1

                    # MAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if magnetizationMode in com:
                        j = magnetizationSteps - 1
                        j_cond = magnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 - (leftReservoir_length +
                                    leftHEXpositions + MCM_length +
                                    rightHEXpositions +
                                    rightReservoir_length) / 2 - steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+n+init_pos), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos+n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos+n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p+1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length/2 - (leftReservoir_length +
                                                      leftHEXpositions +
                                                      MCM_length +
                                                      rightHEXpositions +
                                                      rightReservoir_length) /
                                    2 - steps/2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos+n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1/(2.*freq) -
                                 applied_static_field_time_ratio[1]/2)
                        if n*time_step < applied_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                md = magnetizationMode
                                ap1 = applied_static_field_time_ratio[0]
                                ap2 = applied_static_field_time_ratio[1]
                                ms = magnetizationSteps
                                time_interval = fields.operating_mode(md, ap1,
                                                                      ap2, ms,
                                                                      freq, j)
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val and flag:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)

                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if magnetizationMode not in com:
                                        if j < magnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    # DEMAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if demagnetizationMode in com:
                        j = demagnetizationSteps - 1
                        j_cond = demagnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+init_pos-n), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos-n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] /
                                             (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos-n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p + 1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos-n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1 / (2. * freq) -
                                 removed_static_field_time_ratio[1] / 2)
                        if n*time_step < removed_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                v1 = demagnetizationMode
                                v2 = removed_static_field_time_ratio[0]
                                v3 = removed_static_field_time_ratio[1]
                                v4 = demagnetizationSteps
                                val = fields.operating_mode(v1, v2, v3,
                                                            v4, freq,
                                                            j)
                                time_interval = val
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = demagnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step - previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    ds = demagnetizationSteps
                                    if time_passed == 0. and j < ds and flag:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if demagnetizationMode not in com:
                                        if j < demagnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    if (AMR.q1-q1) == value1 or value1 == 0.:
                        condition = cycle_number < maxCycleNumber - 1
                    else:
                        cond1 = cycle_number < minCycleNumber
                        val = (abs(AMR.objects[2].Q0[leftReservoir_length/2] -
                               value1))
                        cond2 = ((AMR.q1-q1)-value1)/value1 > stopCriteria
                        cond3 = cycle_number < maxCycleNumber
                        condition = (cond1 or cond2) and cond3
                    cycle_number = cycle_number + 1

            elif startingField == 'demagnetization':
                condition = True
                while condition:

                    value1 = AMR.q1 - q1
                    q1 = AMR.q1
                    q2 = AMR.q2

                    if mod_freq != 'default':
                        val = AMR.objects[1].temperature[mod_freq[1]][0]
                        temperature_sensor = val
                        input = open(mod_freq[0], 'r')
                        s = input.readlines()
                        xMod = []
                        yMod = []
                        for line in s:
                            pair = line.split(',')
                            xMod.append(float(pair[0]))
                            yMod.append(float(pair[1]))
                        input.close()
                        freq = np.interp(temperature_sensor, xMod, yMod)
                        steps = int(velocity / (2 * freq * dx))
                        time_step = (1 / (2. * freq)) / steps
                        if int(time_step / dt) == 0:
                            print 'dt or frequency too low'
                        if write_interval > int(time_step / dt):
                            write_interval = 1

                    # DEMAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if demagnetizationMode in com:
                        j = demagnetizationSteps - 1
                        j_cond = demagnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+init_pos-n), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos-n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 + steps/2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] /
                                             (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos-n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p + 1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length / 2 -
                                    (leftReservoir_length + leftHEXpositions +
                                     MCM_length + rightHEXpositions +
                                     rightReservoir_length) / 2 + steps / 2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos-n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1 / (2. * freq) -
                                 removed_static_field_time_ratio[1] / 2)
                        if n*time_step < removed_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                v1 = demagnetizationMode
                                v2 = removed_static_field_time_ratio[0]
                                v3 = removed_static_field_time_ratio[1]
                                v4 = demagnetizationSteps
                                val = fields.operating_mode(v1, v2, v3,
                                                            v4, freq,
                                                            j)
                                time_interval = val
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = demagnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step - previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    ds = demagnetizationSteps
                                    if time_passed == 0. and j < ds and flag:
                                        v1 = demagnetizationMode
                                        v2 = removed_static_field_time_ratio[0]
                                        v3 = removed_static_field_time_ratio[1]
                                        v4 = demagnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 demagnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  demagnetizationSteps + 1)
                                        AMR.objects[1].deactivate(first,
                                                                  second)
                                        j = j + j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if demagnetizationMode not in com:
                                        if j < demagnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    # MAGNETIZATION

                    j = 0
                    time_passed = 0.
                    j_cond = 0
                    j_inc = 1
                    com = ['constant_left', 'accelerated_left',
                           'decelerated_left']
                    if magnetizationMode in com:
                        j = magnetizationSteps - 1
                        j_cond = magnetizationSteps - 1
                        j_inc = -1

                    for n in range(steps):

                        # contacts
                        contacts = set()

                        init_pos = (fluid_length / 2 - (leftReservoir_length +
                                    leftHEXpositions + MCM_length +
                                    rightHEXpositions +
                                    rightReservoir_length) / 2 - steps / 2)
                        for i in range(leftReservoir_length):
                            contacts.add(((0, i+n+init_pos), (2, i+1),
                                          h_leftreservoir_fluid))

                        if mcm_discontinuity == 'default':
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                contacts.add(((0, i+init_pos+n), (1, i+1),
                                              h_mcm_fluid))
                        else:
                            p = 1
                            init_pos = (leftReservoir_length +
                                        leftHEXpositions + fluid_length / 2 -
                                        (leftReservoir_length +
                                         leftHEXpositions + MCM_length +
                                         rightHEXpositions +
                                         rightReservoir_length) / 2 -
                                        steps / 2)
                            for i in range(MCM_length):
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                cond2 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) -
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if not (i+1 < cond1 and i+1 >= cond2):
                                    contacts.add(((0, i+init_pos+n), (1, i+1),
                                                  h_mcm_fluid))
                                cond1 = (p * len(AMR.objects[1].temperature) /
                                         (mcm_discontinuity[0] + 1) +
                                         int(mcm_discontinuity[1] / (2 * dx)))
                                if i+1 == cond1 and p < mcm_discontinuity[0]:
                                    p = p+1

                        init_pos = (leftReservoir_length + leftHEXpositions +
                                    rightHEXpositions + MCM_length +
                                    fluid_length/2 - (leftReservoir_length +
                                                      leftHEXpositions +
                                                      MCM_length +
                                                      rightHEXpositions +
                                                      rightReservoir_length) /
                                    2 - steps/2)
                        for i in range(rightReservoir_length):
                            contacts.add(((0, i+init_pos+n), (3, i+1),
                                          h_rightreservoir_fluid))

                        AMR.contacts = contacts

                        cond1 = (1/(2.*freq) -
                                 applied_static_field_time_ratio[1]/2)
                        if n*time_step < applied_static_field_time_ratio[0]/2:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)

                        elif n*time_step > cond1:
                            AMR.compute(time_step, write_interval,
                                        solver=solverMode)
                        else:

                            # MCE

                            cond = True
                            previous_time = 0.
                            flag = True
                            if j == j_cond:
                                md = magnetizationMode
                                ap1 = applied_static_field_time_ratio[0]
                                ap2 = applied_static_field_time_ratio[1]
                                ms = magnetizationSteps
                                time_interval = fields.operating_mode(md, ap1,
                                                                      ap2, ms,
                                                                      freq, j)
                            while cond:
                                if (time_interval-time_passed) > time_step:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                    else:
                                        delta_t = time_step
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)
                                    time_passed = time_passed + delta_t
                                    cond = False
                                    flag = True
                                else:
                                    val = magnetizationSteps
                                    if time_passed == 0. and j < val and flag:
                                        v1 = magnetizationMode
                                        v2 = applied_static_field_time_ratio[0]
                                        v3 = applied_static_field_time_ratio[1]
                                        v4 = magnetizationSteps
                                        val = fields.operating_mode(v1, v2, v3,
                                                                    v4, freq,
                                                                    j)
                                        time_interval = val
                                        first = (j * MCM_length /
                                                 magnetizationSteps + 1)
                                        second = ((j+1) * MCM_length /
                                                  magnetizationSteps + 1)
                                        AMR.objects[1].activate(first, second)
                                        j = j+j_inc
                                        delta_t = time_step-previous_time
                                        flag = False
                                        time_passed = delta_t
                                    else:
                                        delta_t = time_interval-time_passed
                                        time_passed = 0.
                                        flag = True
                                    previous_time = delta_t
                                    AMR.compute(delta_t, write_interval,
                                                solver=solverMode)

                                    com = ['constant_left', 'accelerated_left',
                                           'decelerated_left']
                                    if magnetizationMode not in com:
                                        if j < magnetizationSteps:
                                            cond = True
                                        else:
                                            cond = False
                                    else:
                                        if j >= 0:
                                            cond = True
                                        else:
                                            cond = False

                    if (AMR.q1-q1) == value1 or value1 == 0.:
                        condition = cycle_number < maxCycleNumber - 1
                    else:
                        cond1 = cycle_number < minCycleNumber
                        val = abs(AMR.objects[2].Q0[leftReservoir_length/2] -
                                  value1)
                        cond2 = ((AMR.q1-q1)-value1)/value1 > stopCriteria
                        cond3 = cycle_number < maxCycleNumber
                        condition = (cond1 or cond2) and cond3
                    cycle_number = cycle_number + 1

            else:
                print 'incorrect startingField'

        cond1 = startingField == 'magnetization'
        cond2 = startingField == 'demagnetization'
        if cond1 or cond2:

            endTime = time.time()
            simulationTime = endTime - startTime
            hours = int(simulationTime / 3600)
            minutes = int((simulationTime - hours * 3600) / 60)
            seconds = int(simulationTime - hours * 3600 - (minutes * 60))
            hours = '%02d' % hours
            minutes = '%02d' % minutes
            seconds = '%02d' % seconds

            print '------------------------------------------------------'
            print ''
            print 'Number of cycles:', cycle_number
            val = (abs(AMR.objects[2].Q0[leftReservoir_length/2]-value1) /
                   value1)
            print 'Final cycle error:', val
            if mode == 'refrigerator':
                if type_study == 'no load':
                    v1 = leftReservoir_length/2
                    v2 = rightReservoir_length/2
                    val = (abs(AMR.objects[2].temperature[v1][0] -
                           AMR.objects[3].temperature[v2][0]))
                    temperature_span = val
                    print 'No load temperature span (K):', temperature_span
                if type_study == 'fixed temperature span':
                    cooling_power = (-AMR.q2+q2)*freq
                    heating_power = (AMR.q1-q1)*freq
                    working_power = (AMR.q2-q2+AMR.q1-q1)*freq
                    COP = cooling_power/working_power
                    print 'Cooling power (W):', cooling_power
                    print 'Heating power (W):', heating_power
                    print 'Working power (W)', working_power
                    print 'COP:', COP
            if mode == 'heat_pump':
                if type_study == 'no load':
                    v1 = leftReservoir_length/2
                    v2 = rightReservoir_length/2
                    val = (abs(AMR.objects[2].temperature[v1][0] -
                           AMR.objects[3].temperature[v2][0]))
                    temperature_span = val
                    print 'No load temperature span (K):', temperature_span
                if type_study == 'fixed temperature span':
                    cooling_power = (-AMR.q2+q2)*freq
                    heating_power = (AMR.q1-q1)*freq
                    working_power = (AMR.q2-q2+AMR.q1-q1)*freq
                    print 'Cooling power (W):', cooling_power
                    print 'Heating power (W):', heating_power
                    print 'Working power (W)', working_power
                    print 'COP:', heating_power/working_power
            print 'Final time (s):', AMR.objects[0].timePassed
            print 'Simulation duration:', hours + ':' + minutes + ':' + seconds
            print ''
            print '------------------------------------------------------'
            print ''
            print ''
            print ''

        else:
            print 'simulation not complete'
