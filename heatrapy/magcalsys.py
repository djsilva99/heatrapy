# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class magcalsys_solidstate_1D.

Used to compute 1-dimensional models for
magnetocaloric systems

"""

import heatcond
import time
import numpy as np


class magcalsys_solidstate_1D:

    """magcalsys_solidstate_1D class

    computes the thermal processes for 1-dimensional fully solid state
    magnetocaloric systems. The active magnetcic regenerative processes
    can be used with the several allowed demagnetization processes.
    Cascades of materials can also be simulated.

    """

    def __init__(self, fileName, ambTemperature=293,
                 leftThermalSwitch_length=10, rightThermalSwitch_length=10,
                 MCM_length=50, rightReservoir_length=15,
                 leftReservoir_length=15, MCM_material='Gd',
                 leftThermalSwitch_material='idealTS_hot',
                 rightThermalSwitch_material='idealTS_cold',
                 leftReservoir_material='Cu', rightReservoir_material='Cu',
                 freq=.1, dt=0.01, dx=0.002, stopCriteria=5e-8,
                 solverMode='implicit_k(x)', minCycleNumber=50,
                 maxCycleNumber=10000, demagnetizationSteps=1,
                 magnetizationSteps=1, demagnetizationMode='constant_right',
                 magnetizationMode='constant_left', cyclePoints=25,
                 boundaries=[0, 0], note=None, temperatureSensor='default',
                 heatPoints='default', mode='refrigerator', version=None,
                 restingTimeHot='default', restingTimeCold='default',
                 startingField='magnetization'):
        """loads the data and computes.

        fileName: file name where the temperature and heat flux are saved
        ambTemperature: ambient temperature of the whole system
        leftThermalSwitch_length: length of the left thermal switch
        rightThermalSwitch_length: length of the right thermal switch
        leftThermalSwitch_material: string for the material of the
            left thermal switch
        rightThermalSwitch_material: string for the material of the
            right thermal switch
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
        boundaries: list with the boundary conditions
        temperatureSensor: list of two space indexes used to determine the
            temperature span at the end of the simulation. The first term is
            the sensor at the hot end and the second at the cold end
        heatPoints: list of two space indexes used to determine the heat
            flux for the hot end (first term) and cold end (second term)
        mode: mode used for the power calculations (e.g. COP) performed at
            the end of the simulation

        """

        # restingTimes definition
        if restingTimeHot == 'default':
            restingTimeHot = 0.

        if restingTimeCold == 'default':
            restingTimeCold = 0.

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
        print 'Module: magcalsys_solidstate_1D'
        if note is not None:
            print ''
            print 'Note:', note
        print ''
        print 'Mode:', mode
        print 'System:', leftReservoir_material + '/' + \
            leftThermalSwitch_material + '/' + MCM_material + '/' + \
            rightThermalSwitch_material + '/' + rightReservoir_material
        print 'Dimensions (m):', str(dx * leftReservoir_length) + '/' + \
            str(dx * leftThermalSwitch_length) + '/' + str(dx * MCM_length) + \
            '/' + str(dx * rightThermalSwitch_length) + '/' + \
            str(dx * rightReservoir_length)
        print 'Number of points:', str(leftReservoir_length) + '/' + \
            str(leftThermalSwitch_length) + '/' + str(MCM_length) + '/' + \
            str(rightThermalSwitch_length) + '/' + str(rightReservoir_length)
        print 'dx (m):', dx
        print 'dt (s):', dt
        print 'Frequency (Hz):', freq
        print 'Hot resting time ratio:', restingTimeHot
        print 'Cold resting time ratio:', restingTimeCold
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

        # used to calculate the simulation time at the end
        startTime = time.time()

        if heatPoints == 'default':
            leftHeatSensor = leftReservoir_length - 3
            rightHeatSensor = -rightReservoir_length + 3
        else:
            leftHeatSensor = heatPoints[0]
            rightHeatSensor = heatPoints[1]

        if startingField != 'magnetization':
            initialState = True
        else:
            initialState = False

        # initializes the object for the simulation

        materials = [leftReservoir_material, leftThermalSwitch_material,
                     MCM_material, rightThermalSwitch_material,
                     rightReservoir_material]

        lRl = leftReservoir_length
        rRl = rightReservoir_length
        ltsl = leftThermalSwitch_length
        rtsl = rightThermalSwitch_length
        borders = [1, lRl + 1,
                   lRl + ltsl + 1,
                   lRl + ltsl + MCM_length + 1,
                   lRl + ltsl + MCM_length + rtsl + 1,
                   lRl + ltsl + MCM_length + rtsl + rRl + 1]
        heatPoints = [leftHeatSensor, rightHeatSensor]

        a = heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                           fileName=fileName,
                                           materials=materials,
                                           borders=borders,
                                           materialsOrder=[0, 1, 2, 3, 4],
                                           boundaries=boundaries,
                                           heatPoints=heatPoints,
                                           initialState=initialState)

        # defines some variable for the cycles
        if temperatureSensor == 'default':
            rightTemperatureSensor = -(rightReservoir_length / 2)
            leftTemperatureSensor = (leftReservoir_length / 2)
        else:
            rightTemperatureSensor = temperatureSensor[1]
            leftTemperatureSensor = temperatureSensor[0]

        value1 = ambTemperature
        value2 = ambTemperature
        i = 0
        period = (1. / freq)
        stopCriteria2 = 0.
        maximumPower = 0.
        maximumWorkingPower = 0.

        # cycle simulation
        while ((abs((value1 - value2) / value2) > stopCriteria or
                i < minCycleNumber or
                abs((value1 - value2) / value2) > stopCriteria2) and
               i < maxCycleNumber):

            stopCriteria2 = abs((value1 - value2) / value2)
            heatLeft = a.heatLeft
            heatRight = a.heatRight

            if startingField == 'magnetization':

                # MAGNETIZATION AND COMPUTATION

                # mode 1
                if magnetizationMode == 'constant_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                time_interval = ((1 - restingTimeHot -
                                                  restingTimeCold) *
                                                 period / 2.)
                                time_interval = (time_interval /
                                                 (magnetizationSteps - 1))
                                a.compute(time_interval,
                                          int(1 / (freq * dt * cyclePoints)),
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            time_interval = ((1 - restingTimeHot -
                                              restingTimeCold) *
                                             period / 2.)
                            time_interval = (time_interval /
                                             (magnetizationSteps))
                            a.compute(time_interval,
                                      int(1 / (freq * dt * cyclePoints)),
                                      solver=solverMode)

                # mode 2
                if magnetizationMode == 'accelerated_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(j + 1) -
                                            np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            totalmsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 3
                if magnetizationMode == 'decelerated_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(magnetizationSteps - j) -
                                            np.sqrt(magnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            totalmsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(magnetizationSteps - j) -
                                        np.sqrt(magnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 4
                if magnetizationMode == 'constant_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                time_interval = (((1 - restingTimeHot -
                                                  restingTimeCold) * period /
                                                 2.) / (magnetizationSteps -
                                                        1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                              restingTimeCold) * period /
                                             2.) / magnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 5
                if magnetizationMode == 'accelerated_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(magnetizationSteps - j) -
                                            np.sqrt(magnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalmsteps = magnetizationSteps - 1
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(magnetizationSteps - j) -
                                        np.sqrt(magnetizationSteps - j - 1)))
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 6
                if magnetizationMode == 'decelerated_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(j + 1) -
                                            np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalmsteps = magnetizationSteps - 1
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                if restingTimeHot != 0.:
                    a.compute(restingTimeHot * period,
                              int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

                # DEMAGNETIZATION AND COMPUTATION

                # mode 1
                if demagnetizationMode == 'constant_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                time_interval = (((1 - restingTimeHot -
                                                   restingTimeCold) *
                                                  period / 2.) /
                                                 (demagnetizationSteps - 1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                               restingTimeCold) *
                                              period / 2.) /
                                             demagnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 2
                if demagnetizationMode == 'accelerated_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(j + 1) - np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalsteps = demagnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. -
                                                        restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 3
                if demagnetizationMode == 'decelerated_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = demagnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(demagnetizationSteps - j) -
                                            np.sqrt(demagnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = demagnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(demagnetizationSteps - j) -
                                        np.sqrt(demagnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 4
                if demagnetizationMode == 'constant_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                time_interval = (((1 - restingTimeHot -
                                                   restingTimeCold) *
                                                  period / 2.) /
                                                 (demagnetizationSteps - 1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                               restingTimeCold) *
                                              period / 2.) /
                                             demagnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 5
                if demagnetizationMode == 'accelerated_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(demagnetizationSteps - j) -
                                            np.sqrt(demagnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = demagnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(demagnetizationSteps - j) -
                                        np.sqrt(demagnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 6
                if demagnetizationMode == 'decelerated_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(j + 1) - np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. -
                                                        restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                if restingTimeCold != 0.:
                    time_interval = restingTimeCold * period
                    write_interval = int(1 / (freq * dt * cyclePoints))
                    a.compute(time_interval, write_interval, solver=solverMode)

            if startingField == 'demagnetization':

                # DEMAGNETIZATION AND COMPUTATION

                # mode 1
                if demagnetizationMode == 'constant_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                time_interval = (((1 - restingTimeHot -
                                                   restingTimeCold) *
                                                  period / 2.) /
                                                 (demagnetizationSteps - 1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                               restingTimeCold) *
                                              period / 2.) /
                                             demagnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 2
                if demagnetizationMode == 'accelerated_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(j + 1) - np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. -
                                                        restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 3
                if demagnetizationMode == 'decelerated_right':
                    for j in range(demagnetizationSteps):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(demagnetizationSteps - j) -
                                            np.sqrt(demagnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = demagnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(demagnetizationSteps - j) -
                                        np.sqrt(demagnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 4
                if demagnetizationMode == 'constant_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                time_interval = (((1 - restingTimeHot -
                                                   restingTimeCold) *
                                                  period / 2.) /
                                                 (demagnetizationSteps - 1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                               restingTimeCold) *
                                              period / 2.) /
                                             demagnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 5
                if demagnetizationMode == 'accelerated_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(demagnetizationSteps - j) -
                                            np.sqrt(demagnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = demagnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(demagnetizationSteps - j) -
                                        np.sqrt(demagnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 6
                if demagnetizationMode == 'decelerated_left':
                    for j in range(demagnetizationSteps - 1, -1, -1):
                        if restingTimeCold != 0.:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            if j != demagnetizationSteps - 1:
                                totalsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq * np.sqrt(totalsteps))) *
                                           (np.sqrt(j + 1) - np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     ((leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      demagnetizationSteps) +
                                     1)
                            second = (leftReservoir_length + (j + 1) *
                                      ((leftThermalSwitch_length + MCM_length +
                                        rightThermalSwitch_length) /
                                       demagnetizationSteps) +
                                      1)
                            a.deactivate(first, second)
                            totalsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. -
                                                        restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                if restingTimeCold != 0.:
                    time_interval = restingTimeCold * period
                    write_interval = int(1 / (freq * dt * cyclePoints))
                    a.compute(time_interval, write_interval, solver=solverMode)

            # MAGNETIZATION AND COMPUTATION

                # mode 1
                if magnetizationMode == 'constant_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     (magnetizationSteps) + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                time_interval = (((1 - restingTimeHot -
                                                  restingTimeCold) * period /
                                                 2.) / (magnetizationSteps -
                                                        1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                              restingTimeCold) * period /
                                             2.) / magnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 2
                if magnetizationMode == 'accelerated_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(j + 1) -
                                            np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            totalmsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 3
                if magnetizationMode == 'decelerated_right':
                    for j in range(magnetizationSteps):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(magnetizationSteps - j) -
                                            np.sqrt(magnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            totalmsteps = magnetizationSteps
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(magnetizationSteps - j) -
                                        np.sqrt(magnetizationSteps - j - 1)))
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 4
                if magnetizationMode == 'constant_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                time_interval = (((1 - restingTimeHot -
                                                  restingTimeCold) * period /
                                                 2.) / (magnetizationSteps -
                                                        1))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(time_interval, write_interval,
                                          solver=solverMode)
                        else:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            time_interval = (((1 - restingTimeHot -
                                              restingTimeCold) * period /
                                             2.) / magnetizationSteps)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(time_interval, write_interval,
                                      solver=solverMode)

                # mode 5
                if magnetizationMode == 'accelerated_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(magnetizationSteps - j) -
                                            np.sqrt(magnetizationSteps - j -
                                                    1)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalmsteps = magnetizationSteps - 1
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(magnetizationSteps - j) -
                                        np.sqrt(magnetizationSteps - j - 1)))
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                # mode 6
                if magnetizationMode == 'decelerated_left':
                    for j in range(magnetizationSteps - 1, -1, -1):
                        if restingTimeHot != 0.:
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            if j != 0:
                                totalmsteps = magnetizationSteps - 1
                                delta_t = ((1 / (2 * (1. / (1. -
                                                            restingTimeHot -
                                                            restingTimeCold)) *
                                                 freq *
                                                 np.sqrt(totalmsteps))) *
                                           (np.sqrt(j + 1) -
                                            np.sqrt(j)))
                                write_interval = int(1 / (freq * dt *
                                                          cyclePoints))
                                a.compute(delta_t, write_interval,
                                          solver=solverMode)
                        else:
                            totalmsteps = magnetizationSteps - 1
                            delta_t = ((1 / (2 * (1. / (1. - restingTimeHot -
                                                        restingTimeCold)) *
                                             freq * np.sqrt(totalmsteps))) *
                                       (np.sqrt(j + 1) - np.sqrt(j)))
                            first = (leftReservoir_length + j *
                                     (leftThermalSwitch_length + MCM_length +
                                      rightThermalSwitch_length) /
                                     magnetizationSteps + 1)
                            second = (leftReservoir_length + (j + 1) *
                                      (leftThermalSwitch_length + MCM_length +
                                       rightThermalSwitch_length) /
                                      magnetizationSteps + 1)
                            a.activate(first, second)
                            write_interval = int(1 / (freq * dt *
                                                      cyclePoints))
                            a.compute(delta_t, write_interval,
                                      solver=solverMode)

                if restingTimeHot != 0.:
                    a.compute(restingTimeHot * period,
                              int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

            if mode == 'heat_pump':
                if -(a.heatLeft - heatLeft) / period > maximumPower:
                    maximumPower = -(a.heatLeft - heatLeft) / period

                working_power = -((a.heatLeft - heatLeft) -
                                  (a.heatRight - heatRight)) / period
                if working_power > maximumWorkingPower:
                    maximumWorkingPower = - \
                        ((a.heatLeft - heatLeft) -
                         (a.heatRight - heatRight)) / period

            if mode == 'refrigerator':
                if -(a.heatRight - heatRight) / period > maximumPower:
                    maximumPower = -(a.heatRight - heatRight) / period

                working_power = -((a.heatLeft - heatLeft) -
                                  (a.heatRight - heatRight)) / period
                if working_power > maximumWorkingPower:
                    maximumWorkingPower = - \
                        ((a.heatLeft - heatLeft) -
                         (a.heatRight - heatRight)) / period

            if mode == 'refrigerator-energetics':
                coolingPower = -(a.heatRight - heatRight) / period
                workingPower = -((a.heatLeft - heatLeft) -
                                 (a.heatRight - heatRight)) / period
                COP = coolingPower / workingPower
                print (a.heatLeft - heatLeft), (a.heatRight - heatRight)

            if mode == 'heat_pump-energetics':
                heatingPower = -(a.heatLeft - heatLeft) / period
                workingPower = -((a.heatLeft - heatLeft) -
                                 (a.heatRight - heatRight)) / period
                COP = heatingPower / workingPower

            # updates the error values and prints information in log file
            value1 = value2
            value2 = a.temperature[rightTemperatureSensor][1]
            i = 1 + i

        # calculates the final temperatures
        before_final_cycle_time = a.timePassed - 1. / freq
        time_accumulated = 0.
        fileTemperature = open(fileName)
        lines = fileTemperature.readlines()
        final_temperature_right = 0.
        o = 1
        while float(lines[o].split(',')[0]) < before_final_cycle_time:
            lines.pop(0)
            o = o + 1
        for o in range(1, len(lines)):
            time_before = float(lines[o - 1].split(',')[0])
            time_now = float(lines[o].split(',')[0])
            sensor = rightTemperatureSensor - 2
            final_temperature_right = ((time_now - time_before) *
                                       float(lines[o].split(',')[sensor]) +
                                       final_temperature_right)
            time_accumulated = time_accumulated + (time_now - time_before)
        final_temperature_right = final_temperature_right / time_accumulated
        fileTemperature.close()

        before_final_cycle_time = a.timePassed - 1. / freq
        time_accumulated = 0.
        fileTemperature = open(fileName)
        lines = fileTemperature.readlines()
        final_temperature_left = 0.
        o = 1
        while float(lines[o].split(',')[0]) < before_final_cycle_time:
            lines.pop(0)
            o = o + 1
        for o in range(1, len(lines)):
            time_before = float(lines[o - 1].split(',')[0])
            time_now = float(lines[o].split(',')[0])
            sensor = leftTemperatureSensor + 1
            final_temperature_left = ((time_now - time_before) *
                                      float(lines[o].split(',')[sensor]) +
                                      final_temperature_left)
            time_accumulated = time_accumulated + (time_now - time_before)
        final_temperature_left = final_temperature_left / time_accumulated
        fileTemperature.close()

        # prints more information in log file
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
        print 'Number of cycles:', i
        print 'Final cycle error:', abs((value1 - value2) / value2)
        if mode == 'refrigerator':
            print 'Maximum cooling power (W):', maximumPower
            print 'Maximum working power (W)', maximumWorkingPower
            temperature_span = (-final_temperature_right +
                                final_temperature_left)
            print 'No load temperature span (K):', temperature_span
        if mode == 'heat_pump':
            print 'Maximum heating power (W):', maximumPower
            print 'Maximum working power (W)', maximumWorkingPower
            temperature_span = final_temperature_right - final_temperature_left
            print 'No load temperature span (K):', temperature_span
        if mode == 'refrigerator-energetics':
            print 'Cooling power (W):', coolingPower
            print 'Working power (W):', workingPower
            print 'COP:', COP
        if mode == 'heat_pump-energetics':
            print 'Heating power (W):', heatingPower
            print 'Working power (W):', workingPower
            print 'COP:', COP
        print 'Final time (s):', a.timePassed
        print 'Simulation duration:', hours + ':' + minutes + ':' + seconds
        print ''
        print '------------------------------------------------------'
        print ''
        print ''
        print ''


class magcalsys_fluidAMR_1D:

    """magcalsys_fluidAMR_1D class

    computes the thermal processes for 1-dimensional AMR systems that uses
    fluids as heat exchanger.

    """

    def __init__(self, fileName, ambTemperature=293, fluid_length=100,
                 MCM_length=20, rightReservoir_length=3,
                 leftReservoir_length=3, MCM_material='Gd',
                 fluid_material='water', leftReservoir_material='Cu',
                 rightReservoir_material='Cu', freq=.5, dt=0.0001, dx=0.004,
                 stopCriteria=5e-8, solverMode='implicit_k(x)',
                 minCycleNumber=50, maxCycleNumber=10000, cyclePoints=25,
                 note=None, temperatureSensor='default', boundaries=[0, 0],
                 mode='refrigerator', version=None, HEXpositions=30,
                 startingField='magnetization', h=5000000.):

        fileNameMCM = fileName[:-4] + '_MCM.txt'
        MCM = heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                             fileName=fileNameMCM,
                                             materials=[MCM_material],
                                             borders=[1, MCM_length + 1],
                                             materialsOrder=[0],
                                             boundaries=[0, 0],
                                             initialState=True)

        fileNameFluid = fileName[:-4] + '_fluid.txt'
        fluid = heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                               fileName=fileNameFluid,
                                               materials=[fluid_material],
                                               borders=[1, fluid_length + 1],
                                               materialsOrder=[0],
                                               boundaries=[0, 0])

        fileNameHHEX = fileName[:-4] + '_HHEX.txt'
        if boundaries[0] != 0:
            materials = [leftReservoir_material]
            borders = [1, leftReservoir_length + 1]
            HHEX = heatcond.heatcond_activemat_1D(boundaries[0], dx=dx, dt=dt,
                                                  fileName=fileNameHHEX,
                                                  materials=materials,
                                                  borders=borders,
                                                  materialsOrder=[0],
                                                  boundaries=[0, 0])
        else:
            materials = [leftReservoir_material]
            borders = [1, leftReservoir_length + 1]
            HHEX = heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                                  fileName=fileNameHHEX,
                                                  materials=materials,
                                                  borders=borders,
                                                  materialsOrder=[0],
                                                  boundaries=[0, 0])

        fileNameCHEX = fileName[:-4] + '_CHEX.txt'
        if boundaries[1] != 0:
            materials = [rightReservoir_material]
            borders = [1, rightReservoir_length + 1]
            CHEX = heatcond.heatcond_activemat_1D(boundaries[1], dx=dx, dt=dt,
                                                  fileName=fileNameCHEX,
                                                  materials=materials,
                                                  borders=borders,
                                                  materialsOrder=[0],
                                                  boundaries=[0, 0])
        else:
            materials = [rightReservoir_material]
            borders = [1, rightReservoir_length + 1]
            CHEX = heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                                  fileName=fileNameCHEX,
                                                  materials=materials,
                                                  borders=borders,
                                                  materialsOrder=[0],
                                                  boundaries=[0, 0])

        for l in range(100):

            MCM.deactivate(1, MCM_length + 1)

            for i in range(2 * HEXpositions - 1, -1, -1):

                for k in range(int(1. / (4 * freq * HEXpositions * dt))):

                    ind = i + fluid_length / 2 - MCM_length / 2 - HEXpositions
                    Q_MCM = [(h * (-MCM.temperature[j][0] +
                                   fluid.temperature[ind + j][0]))
                             for j in range(1, MCM_length + 2)]
                    Q_HHEX = [(h * (-HHEX.temperature[j][0] +
                                    fluid.temperature[i + j][0]))
                              for j in range(1, leftReservoir_length + 2)]
                    ind = (i + fluid_length - 2 * HEXpositions -
                           rightReservoir_length)
                    Q_CHEX = [(h * (-CHEX.temperature[j][0] +
                                    fluid.temperature[ind + j][0]))
                              for j in range(1, rightReservoir_length + 2)]
                    fluid_space = int((fluid_length - (2 * HEXpositions +
                                                       leftReservoir_length +
                                                       rightReservoir_length +
                                                       MCM_length)) / 2.)
                    ind1 = i + fluid_length / 2 - MCM_length / 2 - HEXpositions
                    ind2 = (i + fluid_length - 2 * HEXpositions -
                            rightReservoir_length)
                    Q_fluid = ([0 for j in range(i)] +
                               [h * (HHEX.temperature[j][0] -
                                     fluid.temperature[i + j][0])
                                for j in range(leftReservoir_length)] +
                               [0 for j in range(fluid_space)] +
                               [h * (MCM.temperature[j][0] -
                                     fluid.temperature[ind1 + j][0])
                                for j in range(MCM_length)] +
                               [0 for j in range(fluid_space)] +
                               [h * (CHEX.temperature[j][0] -
                                     fluid.temperature[ind2 + j][0])
                                for j in range(rightReservoir_length)] +
                               [0 for j in range(2 * HEXpositions - i)])

                    Q_MCM = [(Q_MCM[o], o, o + 1) for o in range(len(Q_MCM))]
                    Q_HHEX = [(Q_HHEX[o], o, o + 1)
                              for o in range(len(Q_HHEX))]
                    Q_CHEX = [(Q_CHEX[o], o, o + 1)
                              for o in range(len(Q_CHEX))]
                    Q_fluid = [(Q_fluid[o], o, o + 1)
                               for o in range(len(Q_fluid))]

                    MCM.changeHeatPower(Q0=Q_MCM)
                    HHEX.changeHeatPower(Q0=Q_HHEX)
                    CHEX.changeHeatPower(Q0=Q_CHEX)
                    fluid.changeHeatPower(Q0=Q_fluid)

                    MCM.compute(dt, 1, solver='implicit_k(x)')
                    if boundaries[0] == 0:
                        HHEX.compute(dt, 1, solver='implicit_k(x)')
                    if boundaries[1] == 0:
                        CHEX.compute(dt, 1, solver='implicit_k(x)')
                    fluid.compute(dt, 1, solver='implicit_k(x)')

            MCM.activate(1, MCM_length + 1)

            for i in range(2 * HEXpositions):

                for k in range(int(1. / (4 * freq * HEXpositions * dt))):

                    ind = i + fluid_length / 2 - MCM_length / 2 - HEXpositions
                    Q_MCM = [(h * (-MCM.temperature[j][0] +
                                   fluid.temperature[ind + j][0]))
                             for j in range(1, MCM_length + 2)]
                    Q_HHEX = [(h * (-HHEX.temperature[j][0] +
                                    fluid.temperature[i + j][0]))
                              for j in range(1, leftReservoir_length + 2)]
                    ind = (i + fluid_length - 2 * HEXpositions -
                           rightReservoir_length)
                    Q_CHEX = [(h * (-CHEX.temperature[j][0] +
                                    fluid.temperature[ind + j][0]))
                              for j in range(1, rightReservoir_length + 2)]
                    fluid_space = int((fluid_length - (2 * HEXpositions +
                                                       leftReservoir_length +
                                                       rightReservoir_length +
                                                       MCM_length)) / 2.)
                    ind1 = i + fluid_length / 2 - MCM_length / 2 - HEXpositions
                    ind2 = (i + fluid_length - 2 * HEXpositions -
                            rightReservoir_length)
                    Q_fluid = ([0 for j in range(i)] +
                               [h * (HHEX.temperature[j][0] -
                                     fluid.temperature[i + j][0])
                                for j in range(leftReservoir_length)] +
                               [0 for j in range(fluid_space)] +
                               [h * (MCM.temperature[j][0] -
                                     fluid.temperature[ind1 + j][0])
                                for j in range(MCM_length)] +
                               [0 for j in range(fluid_space)] +
                               [h * (CHEX.temperature[j][0] -
                                     fluid.temperature[ind2 + j][0])
                                for j in range(rightReservoir_length)] +
                               [0 for j in range(2 * HEXpositions - i)])
                    Q_MCM = [(Q_MCM[o], o, o + 1) for o in range(len(Q_MCM))]
                    Q_HHEX = [(Q_HHEX[o], o, o + 1)
                              for o in range(len(Q_HHEX))]
                    Q_CHEX = [(Q_CHEX[o], o, o + 1)
                              for o in range(len(Q_CHEX))]
                    Q_fluid = [(Q_fluid[o], o, o + 1)
                               for o in range(len(Q_fluid))]

                    MCM.changeHeatPower(Q0=Q_MCM)
                    HHEX.changeHeatPower(Q0=Q_HHEX)
                    CHEX.changeHeatPower(Q0=Q_CHEX)
                    fluid.changeHeatPower(Q0=Q_fluid)

                    MCM.compute(dt, 1, solver='implicit_k(x)')
                    if boundaries[0] == 0:
                        HHEX.compute(dt, 1, solver='implicit_k(x)')
                    if boundaries[1] == 0:
                        CHEX.compute(dt, 1, solver='implicit_k(x)')
                    fluid.compute(dt, 1, solver='implicit_k(x)')
