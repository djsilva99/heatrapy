# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class magcalsys_solidstate_1D.

Used to compute 1-dimensional models for solid state systems involving only one
thermal object and caloric effects.

"""

from .. import objects
from .. import fields
import time
import numpy as np


class magcalsys_solidstate_1D:

    """magcalsys_solidstate_1D class

    computes the thermal processes for 1-dimensional fully solid state
    ferroic-based systems. The active regenerative processes can be used with
    the several allowed modes for application and removal of fields. Cascades
    of materials can also be computed.

    """

    def __init__(self, fileName, ambTemperature=293,
                 leftThermalSwitch_length=2, rightThermalSwitch_length=2,
                 MCM_length=20, rightReservoir_length=3,
                 leftReservoir_length=3, MCM_material=((0.002, 'Gd'),),
                 leftThermalSwitch_material='idealTS_hot',
                 rightThermalSwitch_material='idealTS_cold',
                 leftReservoir_material='Cu', rightReservoir_material='Cu',
                 freq=.1, dt=.01, dx=0.002, stopCriteria=5e-3,
                 solverMode='implicit_k(x)', minCycleNumber=30,
                 maxCycleNumber=31, demagnetizationSteps=3,
                 magnetizationSteps=1, demagnetizationMode='accelerated_right',
                 magnetizationMode='accelerated_left', cyclePoints=25,
                 boundaries=[293, 293], note=None, temperatureSensor='default',
                 heatPoints='default', mode='refrigerator', version=None,
                 restingTimeHot='default', restingTimeCold='default',
                 startingField='magnetization',
                 type_study='fixed temperature span', h_left=50000.,
                 h_right=50000.,
                 mod_freq='default'):
        """Initialization.

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
        version: heatrapy version (default is None)
        type_study: 'no load' or 'fixed temperature span'
        h_left: left heat transfer coefficient
        h_right: right heat transfer coefficient
        mod_freq: if not 'default', allows to modulate the frequency according
            to a specific temperature. tuple with first element as filename,
            and second the sensor point.

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
            leftThermalSwitch_material + '/MCM material/' + \
            rightThermalSwitch_material + '/' + rightReservoir_material
        print 'Dimensions (m):', str(dx * leftReservoir_length) + '/' + \
            str(dx * leftThermalSwitch_length) + '/' + str(dx * MCM_length) + \
            '/' + str(dx * rightThermalSwitch_length) + '/' + \
            str(dx * rightReservoir_length)
        print 'Number of points:', str(leftReservoir_length) + '/' + \
            str(leftThermalSwitch_length) + '/' + str(MCM_length) + '/' + \
            str(rightThermalSwitch_length) + '/' + str(rightReservoir_length)
        print 'MCM material:', MCM_material
        print 'Left heat transfer coefficient (W /(m2 K)):', h_left
        print 'Right heat transfer coefficient (W /(m2 K)):', h_right
        print 'dx (m):', dx
        print 'dt (s):', dt
        print 'Frequency (Hz):', freq
        print 'Hot resting time ratio:', restingTimeHot
        print 'Cold resting time ratio:', restingTimeCold
        print 'Solver:', solverMode
        print 'Frequency modulation:', mod_freq
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
            if leftReservoir_length - 3 > 0:
                leftHeatSensor = leftReservoir_length - 3
            else:
                leftHeatSensor = 2
            if -rightReservoir_length + 3 < 0:
                rightHeatSensor = -rightReservoir_length + 3
            else:
                rightHeatSensor = -2
        else:
            leftHeatSensor = heatPoints[0]
            rightHeatSensor = heatPoints[1]

        if startingField != 'magnetization':
            initialState = True
        else:
            initialState = False

        # initializes the object for the simulation
        materials = [leftReservoir_material, leftThermalSwitch_material,
                     MCM_material[0][1], rightThermalSwitch_material,
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

        a = objects.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt,
                                          fileName=fileName,
                                          materials=materials,
                                          borders=borders,
                                          materialsOrder=[0, 1, 2, 3, 4],
                                          boundaries=boundaries,
                                          heatPoints=heatPoints,
                                          initialState=initialState,
                                          h_left=50000., h_right=50000.)

        # defines the material cascade
        k = leftReservoir_length + leftThermalSwitch_length
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
            a.materials.append(mats.calmatpro(tadi, tadd, cpa, cp0, k0, ka,
                                              rho0, rhoa))
            for j in range(k, k+int(MCM_material[i][0]/dx+1)):
                a.materialsIndex[j] = len(a.materials)-1
            k = k+int(MCM_material[i][0]/dx)

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

            if mod_freq != 'default':
                temperature_sensor = a.temperature[mod_freq[1]][0]
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
                period = (1. / freq)

            if startingField == 'magnetization':

                # MAGNETIZATION AND COMPUTATION

                com = ['constant_left', 'accelerated_left', 'decelerated_left']
                if magnetizationMode in com:
                    step_range = range(magnetizationSteps - 1, -1, -1)
                else:
                    step_range = range(magnetizationSteps)

                for j in step_range:

                    first = (leftReservoir_length + j *
                             (leftThermalSwitch_length + MCM_length +
                              rightThermalSwitch_length) /
                             magnetizationSteps + 1)
                    second = (leftReservoir_length + (j + 1) *
                              (leftThermalSwitch_length + MCM_length +
                               rightThermalSwitch_length) /
                              magnetizationSteps + 1)
                    a.activate(first, second)
                    delta_t = fields.operating_mode(magnetizationMode,
                                                    restingTimeCold,
                                                    restingTimeHot,
                                                    magnetizationSteps, freq,
                                                    j)
                    a.compute(delta_t, int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

                a.compute(restingTimeHot * period,
                          int(1 / (freq * dt * cyclePoints)),
                          solver=solverMode)

                # DEMAGNETIZATION AND COMPUTATION

                com = ['constant_left', 'accelerated_left', 'decelerated_left']
                if magnetizationMode in com:
                    step_range = range(demagnetizationSteps - 1, -1, -1)
                else:
                    step_range = range(demagnetizationSteps)

                for j in step_range:
                    first = (leftReservoir_length + j *
                             (leftThermalSwitch_length + MCM_length +
                              rightThermalSwitch_length) /
                             demagnetizationSteps + 1)
                    second = (leftReservoir_length + (j + 1) *
                              (leftThermalSwitch_length + MCM_length +
                               rightThermalSwitch_length) /
                              demagnetizationSteps + 1)
                    a.deactivate(first, second)
                    delta_t = fields.operating_mode(demagnetizationMode,
                                                    restingTimeCold,
                                                    restingTimeHot,
                                                    demagnetizationSteps, freq,
                                                    j)
                    a.compute(delta_t, int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

                time_interval = restingTimeCold * period
                write_interval = int(1 / (freq * dt * cyclePoints))
                a.compute(time_interval, write_interval, solver=solverMode)

            if startingField == 'demagnetization':

                # DEMAGNETIZATION AND COMPUTATION

                com = ['constant_left', 'accelerated_left', 'decelerated_left']
                if magnetizationMode in com:
                    step_range = range(demagnetizationSteps - 1, -1, -1)
                else:
                    step_range = range(demagnetizationSteps)

                for j in step_range:
                    first = (leftReservoir_length + j *
                             (leftThermalSwitch_length + MCM_length +
                              rightThermalSwitch_length) /
                             demagnetizationSteps + 1)
                    second = (leftReservoir_length + (j + 1) *
                              (leftThermalSwitch_length + MCM_length +
                               rightThermalSwitch_length) /
                              demagnetizationSteps + 1)
                    a.deactivate(first, second)
                    delta_t = fields.operating_mode(demagnetizationMode,
                                                    restingTimeCold,
                                                    restingTimeHot,
                                                    demagnetizationSteps, freq,
                                                    j)
                    a.compute(delta_t, int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

                time_interval = restingTimeCold * period
                write_interval = int(1 / (freq * dt * cyclePoints))
                a.compute(time_interval, write_interval, solver=solverMode)

            # MAGNETIZATION AND COMPUTATION

                com = ['constant_left', 'accelerated_left', 'decelerated_left']
                if magnetizationMode in com:
                    step_range = range(magnetizationSteps - 1, -1, -1)
                else:
                    step_range = range(magnetizationSteps)

                for j in step_range:
                    first = (leftReservoir_length + j *
                             (leftThermalSwitch_length + MCM_length +
                              rightThermalSwitch_length) /
                             magnetizationSteps + 1)
                    second = (leftReservoir_length + (j + 1) *
                              (leftThermalSwitch_length + MCM_length +
                               rightThermalSwitch_length) /
                              magnetizationSteps + 1)
                    a.activate(first, second)
                    delta_t = fields.operating_mode(magnetizationMode,
                                                    restingTimeCold,
                                                    restingTimeHot,
                                                    magnetizationSteps, freq,
                                                    j)
                    a.compute(delta_t, int(1 / (freq * dt * cyclePoints)),
                              solver=solverMode)

                a.compute(restingTimeHot * period,
                          int(1 / (freq * dt * cyclePoints)),
                          solver=solverMode)

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
        if type_study == 'fixed temperature span':
            heating_power = -(a.heatLeft - heatLeft) / period
            cooling_power = -(a.heatRight - heatRight) / period
            working_power = heating_power - cooling_power
            print 'Heating power (W):', heating_power
            print 'Cooling power (W):', cooling_power
            print 'Working power (W):', working_power
            if mode == 'refrigerator':
                print 'COP:', cooling_power / working_power
            if mode == 'heat_pump':
                print 'COP:', heating_power / working_power
        else:
            temperature_span = (-final_temperature_right +
                                final_temperature_left)
            if mode == 'refrigerator':
                print 'No load temperature span (K):', temperature_span
            if mode == 'heat_pump':
                print 'No load temperature span (K):', -temperature_span
        print 'Final time (s):', a.timePassed
        print 'Simulation duration:', hours + ':' + minutes + ':' + seconds
        print ''
        print '------------------------------------------------------'
        print ''
        print ''
        print ''
