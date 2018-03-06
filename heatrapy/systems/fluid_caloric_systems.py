# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the function fluid_active_regenerator.

Used to compute 1-dimensional models for ferroic-based systems using fluids
with the regenerative heat process.

"""

from .. import objects
import time
import numpy as np
from .. import fields


def fluid_active_regenerator(file_name, amb_temperature=298, fluid_length=160,
                             MCM_length=50, right_reservoir_length=20,
                             left_reservoir_length=20,
                             MCM_material=((0.000, 'Gd'),),
                             fluid_material='water',
                             left_reservoir_material='Cu',
                             right_reservoir_material='Cu', freq=0.17, dt=0.1,
                             dx=0.001, stop_criteria=1e-6,
                             solver='implicit_k(x)', min_cycle_number=1,
                             max_cycle_number=5000, cycle_points=25, note=None,
                             boundaries=((2, 298), (3, 0)), mode='heat_pump',
                             version=None, leftHEXpositions=10,
                             rightHEXpositions=10, starting_field='applied',
                             temperature_sensor=(3, -3), field_removal_steps=1,
                             field_applied_steps=1,
                             field_removal_mode='accelerated_left',
                             field_applied_mode='accelerated_right',
                             applied_static_field_time_ratio=(0., 0.),
                             removed_static_field_time_ratio=(0., 0.),
                             h_mcm_fluid=1500000,
                             h_leftreservoir_fluid=1000000,
                             h_rightreservoir_fluid=1000000,
                             mcm_discontinuity='default',
                             type_study='no_load', velocity=.007,
                             mod_freq='default'):
    """fluid_active_regenerator class

    Computes the thermal processes for 1-dimensional ferroic-based systems
    using fluids. The active regenerative processes can be used with the
    several allowed modes for application and removal of fields. Cascades of
    materials can also be computed.

    file_name: file name where the temperature and heat flux are saved
    amb_temperature: ambient temperature of the whole system
    fluid_length: length of the fluid
    fluid_material: string for the material of the fluid
    MCM_length: length of the magnetocaloric material
    MCM_material: string for the material of the magnetocaloric material
    left_reservoir_length: length of the left reservoir
    right_reservoir_length: length of the right reservoir
    left_reservoir_material: string for the material of the left reservoir
    right_reservoir_material: string for the material of the right reservoir
    freq: operating frequency
    dt: times step
    dx: space step
    stop_criteria: error threshold to stop the simulation
    solver: the solver
    min_cycle_number: minimum number of cycles that has to be computed
    max_cycle_number: maximum number of cycles that has to be computed
    field_removal_steps: number of steps during the field removal
    field_applied_steps: number of steps during the application of field
    field_removal_mode: mode of the field removal
        modes can be constant_right, constant_left, accelerated_right,
        accelerated_left, decelerated_right, and decelerated_left
    field_applied_mode is the mode of the application of field
        modes can be constant_right, constant_left, accelerated_right,
        accelerated_left, decelerated_right, and decelerated_left
    cycle_points: number of points recorded for each position for each cycle
    boundaries: tuple of two entries that define the boundary condition
        for tempreture. The first corresponds to the thermal obect while
        the second defines the temperature. If 0 the boundary condition is
        insulation
    temperature_sensor: list of two space indexes used to determine the
        temperature span at the end of the simulation. The first term is
        the sensor at the hot end and the second at the cold end
    mode: mode used for the power calculations (e.g. COP) performed at
        the end of the simulation
    version: heatrapy version (default is None)
    type_study: 'no_load' or 'fixed_temperature_span'
    h_mcm_fluid: heat transfer coefficient for fluid - MCM
    h_leftreservoir_fluid: heat transfer coefficient for
        fluid - left reservoir
    h_rigthreservoir_fluid: heat transfer coefficient for
        fluid - right reservoir
    mod_freq: if not 'default', allows to modulate the frequency according
        to a specific temperature. tuple with first element as file_name,
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

    # check the validity of inputs
    cond01 = isinstance(amb_temperature, float)
    cond01 = cond01 or isinstance(amb_temperature, int)
    cond02 = isinstance(fluid_length, int)
    cond03 = isinstance(velocity, int) or isinstance(velocity, float)
    cond04 = isinstance(MCM_length, int)
    cond05 = isinstance(dx, int) or isinstance(dx, float)
    cond06 = isinstance(dt, int) or isinstance(dt, float)
    cond07 = isinstance(file_name, unicode)
    cond07 = cond07 or isinstance(file_name, str)
    cond08 = isinstance(boundaries, tuple)
    cond09 = isinstance(right_reservoir_length, int)
    cond10 = isinstance(left_reservoir_length, int)
    cond11 = isinstance(MCM_material, tuple)
    cond12 = isinstance(fluid_material, str)
    cond12 = cond12 or isinstance(fluid_material, unicode)
    cond13 = isinstance(left_reservoir_material, str)
    cond13 = cond13 or isinstance(left_reservoir_material, unicode)
    cond14 = isinstance(right_reservoir_material, str)
    cond14 = cond14 or isinstance(right_reservoir_material, unicode)
    cond15 = isinstance(freq, int) or isinstance(freq, float)
    cond16 = isinstance(stop_criteria, int) or isinstance(stop_criteria, float)
    cond17 = isinstance(solver, unicode) or isinstance(solver, str)
    cond18 = isinstance(min_cycle_number, int)
    cond19 = isinstance(max_cycle_number, int)
    cond20 = isinstance(cycle_points, int)
    cond21 = isinstance(mode, unicode) or isinstance(mode, str)
    cond22 = isinstance(leftHEXpositions, int)
    cond23 = isinstance(rightHEXpositions, int)
    cond24 = starting_field == 'applied' or starting_field == 'removal'
    cond25 = isinstance(temperature_sensor, tuple)
    cond26 = isinstance(field_removal_steps, int)
    cond27 = isinstance(field_applied_steps, int)
    allowed_modes = ['constant_right', 'constant_left', 'accelerated_right',
                     'accelerated_left', 'decelerated_right',
                     'decelerated_left']
    cond28 = field_removal_mode in allowed_modes
    cond29 = field_applied_mode in allowed_modes
    cond30 = isinstance(applied_static_field_time_ratio, tuple)
    cond31 = isinstance(removed_static_field_time_ratio, tuple)
    cond32 = isinstance(h_mcm_fluid, int) or isinstance(h_mcm_fluid, float)
    cond33 = isinstance(h_leftreservoir_fluid, int)
    cond33 = cond33 or isinstance(h_leftreservoir_fluid, float)
    cond34 = isinstance(h_rightreservoir_fluid, int)
    cond34 = cond33 or isinstance(h_rightreservoir_fluid, float)
    cond35 = type_study == 'no_load' or type_study == 'fixed_temperature_span'
    condition = cond01 and cond02 and cond03 and cond04 and cond05
    condition = condition and cond06 and cond07 and cond08 and cond09
    condition = condition and cond10 and cond11 and cond12 and cond13
    condition = condition and cond14 and cond15 and cond16 and cond17
    condition = condition and cond18 and cond19 and cond20 and cond21
    condition = condition and cond22 and cond23 and cond24 and cond25
    condition = condition and cond26 and cond27 and cond28 and cond29
    condition = condition and cond30 and cond31 and cond32 and cond33
    condition = condition and cond34 and cond35
    if not condition:
        raise ValueError

    # initial definitions
    if starting_field == 'removal':
        initial_state = True
    else:
        initial_state = False

    cycle_number = 0
    AMR = objects.system_objects(number_objects=4,
                                 materials=(fluid_material,
                                            MCM_material[0][1],
                                            left_reservoir_material,
                                            right_reservoir_material),
                                 objects_length=(fluid_length, MCM_length,
                                                 left_reservoir_length,
                                                 right_reservoir_length),
                                 amb_temperature=amb_temperature, dx=dx, dt=dt,
                                 file_name=file_name,
                                 initial_state=initial_state,
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
            '/../database/vacuum/tadi.txt'
        tadd = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/tadd.txt'
        cpa = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/cpa.txt'
        cp0 = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/cp0.txt'
        k0 = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/k0.txt'
        ka = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/ka.txt'
        rho0 = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/rho0.txt'
        rhoa = os.path.dirname(os.path.realpath(__file__)) + \
            '/../database/vacuum/rhoa.txt'
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

    write_interval = cycle_points/2
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
    print file_name
    print '------------------------------------------------------'
    print ''
    print 'heatconpy version:', version
    print 'Module: fluid_active_regenerator'
    if note is not None:
        print ''
        print 'Note:', note
    print ''
    print 'Mode:', mode
    print 'Fluid:', fluid_material + ' (' + str(fluid_length*dx) + ')'
    print 'MCM material:', MCM_material
    string = ' (' + str(left_reservoir_length * dx) + ')'
    print 'Left reservoir:', left_reservoir_material + string
    string = ' (' + str(right_reservoir_length * dx) + ')'
    print 'Right reservoir:', right_reservoir_material + string
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
    print 'Solver:', solver
    print 'Applied field mode:', field_applied_mode
    print 'Applied field steps:', field_applied_steps
    print 'Field removal mode:', field_removal_mode
    print 'Field removal steps:', field_removal_steps
    print 'Starting Field:', starting_field
    print 'Boundaries:', boundaries
    print 'Ambient temperature (K):', amb_temperature
    print 'Stop criteria:', stop_criteria
    print 'Time:', time.ctime()
    print ''

    start_time = time.time()

    if type_study == 'no_load':

        value1 = AMR.objects[0].temperature[temperature_sensor[0]][0]
        value2 = AMR.objects[0].temperature[temperature_sensor[1]][0]

        if starting_field == 'applied':

            condition = True
            while condition:

                vl1 = left_reservoir_length/2
                val = AMR.objects[2].temperature[vl1][0]
                value1 = val
                vl2 = right_reservoir_length/2
                val = AMR.objects[3].temperature[vl2][0]
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

                # APPLICATION OF FIELD

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_applied_mode in com:
                    j = field_applied_steps - 1
                    j_cond = field_applied_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 - (left_reservoir_length +
                                leftHEXpositions + MCM_length +
                                rightHEXpositions +
                                right_reservoir_length) / 2 - steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+n+init_pos), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
                                    steps / 2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos+n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length/2 - (left_reservoir_length +
                                                  leftHEXpositions +
                                                  MCM_length +
                                                  rightHEXpositions +
                                                  right_reservoir_length) /
                                2 - steps/2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos+n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1/(2.*freq) -
                             applied_static_field_time_ratio[1]/2)
                    if n*time_step < applied_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            md = field_applied_mode
                            ap1 = applied_static_field_time_ratio[0]
                            ap2 = applied_static_field_time_ratio[1]
                            ms = field_applied_steps
                            time_interval = fields.operating_mode(md, ap1,
                                                                  ap2, ms,
                                                                  freq, j)
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_applied_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
                                    AMR.objects[1].activate(first, second)
                                    j = j + j_inc
                                    delta_t = time_step-previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                val = field_applied_steps
                                if time_passed == 0. and j < val and flag:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
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
                                            solver=solver)

                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_applied_mode not in com:
                                    if j < field_applied_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                # FIELD REMOVAL

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_removal_mode in com:
                    j = field_removal_steps - 1
                    j_cond = field_removal_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos-n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1 / (2. * freq) -
                             removed_static_field_time_ratio[1] / 2)
                    if n*time_step < removed_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            v1 = field_removal_mode
                            v2 = removed_static_field_time_ratio[0]
                            v3 = removed_static_field_time_ratio[1]
                            v4 = field_removal_steps
                            val = fields.operating_mode(v1, v2, v3,
                                                        v4, freq,
                                                        j)
                            time_interval = val
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_removal_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
                                    AMR.objects[1].deactivate(first,
                                                              second)
                                    j = j + j_inc
                                    delta_t = time_step - previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                ds = field_removal_steps
                                if time_passed == 0. and j < ds and flag:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
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
                                            solver=solver)
                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_removal_mode not in com:
                                    if j < field_removal_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                if value1 == value2:
                    condition = cycle_number < max_cycle_number - 1
                else:
                    cond1 = cycle_number < min_cycle_number
                    vl1 = left_reservoir_length/2
                    vl2 = right_reservoir_length/2
                    cond2 = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                      AMR.objects[3].temperature[vl2][0]) -
                                  abs(value1-value2)))/abs(value1-value2))
                    cond2 = cond2 > stop_criteria
                    cond3 = cycle_number < max_cycle_number - 1
                    condition = (cond1 or cond2) and cond3

                cycle_number = cycle_number + 1

        elif starting_field == 'removal':
            condition = True
            while condition:
                vl1 = left_reservoir_length/2
                val = AMR.objects[2].temperature[vl1][0]
                value1 = val
                vl2 = right_reservoir_length/2
                val = AMR.objects[3].temperature[vl2][0]
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
                if field_removal_mode in com:
                    j = field_removal_steps - 1
                    j_cond = field_removal_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) /
                                    2 + steps / 2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos-n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 +
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts
                    val = (1 / (2. * freq) -
                           removed_static_field_time_ratio[1] / 2)
                    if n*time_step < removed_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > val:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            v1 = field_removal_mode
                            v2 = removed_static_field_time_ratio[0]
                            v3 = removed_static_field_time_ratio[1]
                            v4 = field_removal_steps
                            val = fields.operating_mode(v1, v2, v3,
                                                        v4, freq,
                                                        j)
                            time_interval = val
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_removal_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
                                    AMR.objects[1].deactivate(first,
                                                              second)
                                    j = j+j_inc
                                    delta_t = time_step - previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                val = field_removal_steps
                                if time_passed == 0. and j < val and flag:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
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
                                            solver=solver)
                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_removal_mode not in com:
                                    if j < field_removal_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                # APPLICATION OF FIELD

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_applied_mode in com:
                    j = field_applied_steps - 1
                    j_cond = field_applied_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 - (left_reservoir_length +
                                leftHEXpositions + MCM_length +
                                rightHEXpositions +
                                right_reservoir_length) / 2 - steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+n+init_pos), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
                                    steps / 2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos+n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length/2 - (left_reservoir_length +
                                                  leftHEXpositions +
                                                  MCM_length +
                                                  rightHEXpositions +
                                                  right_reservoir_length) /
                                2 - steps/2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos+n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1/(2.*freq) -
                             applied_static_field_time_ratio[1]/2)
                    if n*time_step < applied_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            md = field_applied_mode
                            ap1 = applied_static_field_time_ratio[0]
                            ap2 = applied_static_field_time_ratio[1]
                            ms = field_applied_steps
                            time_interval = fields.operating_mode(md, ap1,
                                                                  ap2, ms,
                                                                  freq, j)
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_applied_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
                                    AMR.objects[1].activate(first, second)
                                    j = j+j_inc
                                    delta_t = time_step-previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                val = field_applied_steps
                                if time_passed == 0. and j < val and flag:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
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
                                            solver=solver)

                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_applied_mode not in com:
                                    if j < field_applied_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                if value1 == value2:
                    condition = cycle_number < max_cycle_number - 1
                else:
                    cond1 = cycle_number < min_cycle_number
                    vl1 = left_reservoir_length/2
                    vl2 = right_reservoir_length/2
                    cond2 = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                      AMR.objects[3].temperature[vl2][0]) -
                                  abs(value1-value2)))/abs(value1-value2))
                    cond2 = cond2 > stop_criteria
                    cond3 = cycle_number < max_cycle_number - 1
                    condition = (cond1 or cond2) and cond3

                cycle_number = cycle_number + 1

        else:
            print 'incorrect starting_field'

    if type_study == 'fixed_temperature_span':

        if starting_field == 'applied':

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

                # APPLICATION OF FIELD

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_applied_mode in com:
                    j = field_applied_steps - 1
                    j_cond = field_applied_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 - (left_reservoir_length +
                                leftHEXpositions + MCM_length +
                                rightHEXpositions +
                                right_reservoir_length) / 2 - steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+n+init_pos), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
                                    steps / 2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos+n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length/2 - (left_reservoir_length +
                                                  leftHEXpositions +
                                                  MCM_length +
                                                  rightHEXpositions +
                                                  right_reservoir_length) /
                                2 - steps/2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos+n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1/(2.*freq) -
                             applied_static_field_time_ratio[1]/2)
                    if n*time_step < applied_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            md = field_applied_mode
                            ap1 = applied_static_field_time_ratio[0]
                            ap2 = applied_static_field_time_ratio[1]
                            ms = field_applied_steps
                            time_interval = fields.operating_mode(md, ap1,
                                                                  ap2, ms,
                                                                  freq, j)
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_applied_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
                                    AMR.objects[1].activate(first, second)
                                    j = j+j_inc
                                    delta_t = time_step-previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                val = field_applied_steps
                                if time_passed == 0. and j < val and flag:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
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
                                            solver=solver)

                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_applied_mode not in com:
                                    if j < field_applied_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                # FIELD REMOVAL

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_removal_mode in com:
                    j = field_removal_steps - 1
                    j_cond = field_removal_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos-n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1 / (2. * freq) -
                             removed_static_field_time_ratio[1] / 2)
                    if n*time_step < removed_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            v1 = field_removal_mode
                            v2 = removed_static_field_time_ratio[0]
                            v3 = removed_static_field_time_ratio[1]
                            v4 = field_removal_steps
                            val = fields.operating_mode(v1, v2, v3,
                                                        v4, freq,
                                                        j)
                            time_interval = val
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_removal_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
                                    AMR.objects[1].deactivate(first,
                                                              second)
                                    j = j + j_inc
                                    delta_t = time_step - previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                ds = field_removal_steps
                                if time_passed == 0. and j < ds and flag:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
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
                                            solver=solver)
                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_removal_mode not in com:
                                    if j < field_removal_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                if (AMR.q1-q1) == value1 or value1 == 0.:
                    condition = cycle_number < max_cycle_number - 1
                else:
                    cond1 = cycle_number < min_cycle_number
                    val = (abs(AMR.objects[2].Q0[left_reservoir_length/2] -
                           value1))
                    cond2 = ((AMR.q1-q1)-value1)/value1 > stop_criteria
                    cond3 = cycle_number < max_cycle_number
                    condition = (cond1 or cond2) and cond3
                cycle_number = cycle_number + 1

        elif starting_field == 'removal':
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

                # FIELD REMOVAL

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_removal_mode in com:
                    j = field_removal_steps - 1
                    j_cond = field_removal_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos-n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 + steps/2)
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length / 2 -
                                (left_reservoir_length + leftHEXpositions +
                                    MCM_length + rightHEXpositions +
                                    right_reservoir_length) / 2 + steps / 2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos-n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1 / (2. * freq) -
                             removed_static_field_time_ratio[1] / 2)
                    if n*time_step < removed_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            v1 = field_removal_mode
                            v2 = removed_static_field_time_ratio[0]
                            v3 = removed_static_field_time_ratio[1]
                            v4 = field_removal_steps
                            val = fields.operating_mode(v1, v2, v3,
                                                        v4, freq,
                                                        j)
                            time_interval = val
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_removal_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
                                    AMR.objects[1].deactivate(first,
                                                              second)
                                    j = j + j_inc
                                    delta_t = time_step - previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                ds = field_removal_steps
                                if time_passed == 0. and j < ds and flag:
                                    v1 = field_removal_mode
                                    v2 = removed_static_field_time_ratio[0]
                                    v3 = removed_static_field_time_ratio[1]
                                    v4 = field_removal_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_removal_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_removal_steps + 1)
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
                                            solver=solver)
                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_removal_mode not in com:
                                    if j < field_removal_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                # APPLICATION OF FIELD

                j = 0
                time_passed = 0.
                j_cond = 0
                j_inc = 1
                com = ['constant_left', 'accelerated_left',
                       'decelerated_left']
                if field_applied_mode in com:
                    j = field_applied_steps - 1
                    j_cond = field_applied_steps - 1
                    j_inc = -1

                for n in range(steps):

                    # contacts
                    contacts = set()

                    init_pos = (fluid_length / 2 - (left_reservoir_length +
                                leftHEXpositions + MCM_length +
                                rightHEXpositions +
                                right_reservoir_length) / 2 - steps / 2)
                    for i in range(left_reservoir_length):
                        contacts.add(((0, i+n+init_pos), (2, i+1),
                                      h_leftreservoir_fluid))

                    if mcm_discontinuity == 'default':
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
                                    steps / 2)
                        for i in range(MCM_length):
                            contacts.add(((0, i+init_pos+n), (1, i+1),
                                          h_mcm_fluid))
                    else:
                        p = 1
                        init_pos = (left_reservoir_length +
                                    leftHEXpositions + fluid_length / 2 -
                                    (left_reservoir_length +
                                        leftHEXpositions + MCM_length +
                                        rightHEXpositions +
                                        right_reservoir_length) / 2 -
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

                    init_pos = (left_reservoir_length + leftHEXpositions +
                                rightHEXpositions + MCM_length +
                                fluid_length/2 - (left_reservoir_length +
                                                  leftHEXpositions +
                                                  MCM_length +
                                                  rightHEXpositions +
                                                  right_reservoir_length) /
                                2 - steps/2)
                    for i in range(right_reservoir_length):
                        contacts.add(((0, i+init_pos+n), (3, i+1),
                                      h_rightreservoir_fluid))

                    AMR.contacts = contacts

                    cond1 = (1/(2.*freq) -
                             applied_static_field_time_ratio[1]/2)
                    if n*time_step < applied_static_field_time_ratio[0]/2:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)

                    elif n*time_step > cond1:
                        AMR.compute(time_step, write_interval,
                                    solver=solver)
                    else:

                        # MCE

                        cond = True
                        previous_time = 0.
                        flag = True
                        if j == j_cond:
                            md = field_applied_mode
                            ap1 = applied_static_field_time_ratio[0]
                            ap2 = applied_static_field_time_ratio[1]
                            ms = field_applied_steps
                            time_interval = fields.operating_mode(md, ap1,
                                                                  ap2, ms,
                                                                  freq, j)
                        while cond:
                            if (time_interval-time_passed) > time_step:
                                val = field_applied_steps
                                if time_passed == 0. and j < val:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
                                    AMR.objects[1].activate(first, second)
                                    j = j+j_inc
                                    delta_t = time_step-previous_time
                                else:
                                    delta_t = time_step
                                AMR.compute(delta_t, write_interval,
                                            solver=solver)
                                time_passed = time_passed + delta_t
                                cond = False
                                flag = True
                            else:
                                val = field_applied_steps
                                if time_passed == 0. and j < val and flag:
                                    v1 = field_applied_mode
                                    v2 = applied_static_field_time_ratio[0]
                                    v3 = applied_static_field_time_ratio[1]
                                    v4 = field_applied_steps
                                    val = fields.operating_mode(v1, v2, v3,
                                                                v4, freq,
                                                                j)
                                    time_interval = val
                                    first = (j * MCM_length /
                                             field_applied_steps + 1)
                                    second = ((j+1) * MCM_length /
                                              field_applied_steps + 1)
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
                                            solver=solver)

                                com = ['constant_left', 'accelerated_left',
                                       'decelerated_left']
                                if field_applied_mode not in com:
                                    if j < field_applied_steps:
                                        cond = True
                                    else:
                                        cond = False
                                else:
                                    if j >= 0:
                                        cond = True
                                    else:
                                        cond = False

                if (AMR.q1-q1) == value1 or value1 == 0.:
                    condition = cycle_number < max_cycle_number - 1
                else:
                    cond1 = cycle_number < min_cycle_number
                    val = abs(AMR.objects[2].Q0[left_reservoir_length/2] -
                              value1)
                    cond2 = ((AMR.q1-q1)-value1)/value1 > stop_criteria
                    cond3 = cycle_number < max_cycle_number
                    condition = (cond1 or cond2) and cond3
                cycle_number = cycle_number + 1

        else:
            print 'incorrect starting_field'

    cond1 = starting_field == 'applied'
    cond2 = starting_field == 'removal'
    if cond1 or cond2:

        endTime = time.time()
        simulationTime = endTime - start_time
        hours = int(simulationTime / 3600)
        minutes = int((simulationTime - hours * 3600) / 60)
        seconds = int(simulationTime - hours * 3600 - (minutes * 60))
        hours = '%02d' % hours
        minutes = '%02d' % minutes
        seconds = '%02d' % seconds

        print '------------------------------------------------------'
        print ''
        if mode == 'refrigerator':
            if type_study == 'no_load':
                v1 = left_reservoir_length/2
                v2 = right_reservoir_length/2
                val = (abs(AMR.objects[2].temperature[v1][0] -
                       AMR.objects[3].temperature[v2][0]))
                temperature_span = val
                error = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                  AMR.objects[3].temperature[vl2][0]) -
                              abs(value1-value2)))/abs(value1-value2))
                print 'Final cycle error:', error
                print 'No load temperature span (K):', temperature_span
            if type_study == 'fixed_temperature_span':
                cooling_power = (-AMR.q2+q2)*freq
                heating_power = (AMR.q1-q1)*freq
                working_power = (AMR.q2-q2+AMR.q1-q1)*freq
                COP = cooling_power/working_power
                error = (abs(AMR.objects[2].Q0[left_reservoir_length/2] -
                             value1) / value1)
                print 'Final cycle error:', error
                print 'Cooling power (W):', cooling_power
                print 'Heating power (W):', heating_power
                print 'Working power (W)', working_power
                print 'COP:', COP
        if mode == 'heat_pump':
            if type_study == 'no_load':
                v1 = left_reservoir_length/2
                v2 = right_reservoir_length/2
                val = (abs(AMR.objects[2].temperature[v1][0] -
                       AMR.objects[3].temperature[v2][0]))
                temperature_span = val
                error = ((abs(abs(AMR.objects[2].temperature[vl1][0] -
                                  AMR.objects[3].temperature[vl2][0]) -
                              abs(value1-value2)))/abs(value1-value2))
                print 'Final cycle error:', error
                print 'No load temperature span (K):', temperature_span
            if type_study == 'fixed_temperature_span':
                cooling_power = (-AMR.q2+q2)*freq
                heating_power = (AMR.q1-q1)*freq
                working_power = (AMR.q2-q2+AMR.q1-q1)*freq
                error = (abs(AMR.objects[2].Q0[left_reservoir_length/2] -
                             value1) / value1)
                print 'Final cycle error:', error
                print 'Cooling power (W):', cooling_power
                print 'Heating power (W):', heating_power
                print 'Working power (W)', working_power
                print 'COP:', heating_power/working_power
        print 'Final time (s):', AMR.objects[0].time_passed
        print 'Simulation duration:', hours + ':' + minutes + ':' + seconds
        print ''
        print '------------------------------------------------------'
        print ''
        print ''
        print ''

    else:
        print 'simulation not complete'
