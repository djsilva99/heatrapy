import heatcond
import time
import numpy as np







class magcalsys_solidstate_1D:





	"""
	This class computes the thermal processes for 1-dimensional magnetocaloric systems. 
	The active magnetcic regenerative processes can be used with the several allowed demagnetization processes. 
	Cascades of materials can also be simulated.
	"""





	def __init__(self, fileName, ambTemperature=293, leftThermalSwitch_length=10, rightThermalSwitch_length=10,
		MCM_length=50, rightReservoir_length=15, leftReservoir_length=15, MCM_material='Gd',
		leftThermalSwitch_material='idealTS_hot', rightThermalSwitch_material='idealTS_cold', leftReservoir_material='Cu',
		rightReservoir_material='Cu', freq=.1, dt=0.01, dx=0.002, stopCriteria=5e-8, solverMode='implicit_k(x)',
		minCycleNumber=50, maxCycleNumber=10000, demagnetizationSteps=1, magnetizationSteps=1,
		demagnetizationMode='constant_right', magnetizationMode='constant_left', cyclePoints=25, boundaries=[0,0],
		note=None, temperatureSensor='default', heatPoints='default', mode='refrigerator', version=None, 
		restingTimeHot='default', restingTimeCold='default', startingField='magnetization'):



		"""
		fileName is the file name where the temperature and heat flux are saved
		ambTemperature is the ambient temperature of the whole system
		leftThermalSwitch_length is the length of the left thermal switch
		rightThermalSwitch_length is the length of the right thermal switch
		leftThermalSwitch_material is the string for the material of the left thermal switch
		rightThermalSwitch_material is the string for the material of the right thermal switch
		MCM_length is the length of the magnetocaloric material
		MCM_material is the string for the material of the magnetocaloric material
		leftReservoir_length is the length of the left reservoir
		rightReservoir_length is the length of the right reservoir
		leftReservoir_material is the string for the material of the left reservoir
		rightReservoir_material is the string for the material of the right reservoir
		freq is the operating frequency
		dt is the times step
		dx is the space step
		stopCriteria is the error threshold to stop the simulation
		solverMode is the solver
		minCycleNumber is the minimum number of cycles that has to be computed
		maxCycleNumber is the maximum number of cycles that has to be computed
		demagnetizationSteps is the number of steps during the demagnetization
		magnetizationSteps is the number of steps during the magnetization
		demagnetizationMode is the mode of demagnetization
			modes can be constant_right, constant_left, accelerated_right, accelerated_left, decelerated_right, 
			and decelerated_left
		magnetizationMode is the mode of magnetization
			modes can be constant_right, constant_left, accelerated_right, accelerated_left, decelerated_right, 
			and decelerated_left
		cyclePoints is the number of points recorded for each position for each cycle
		boundaries is the list with the boundary conditions
		temperatureSensor is a list of two space indexes used to determine the temperature span at the end of the 
			simulation. The first term is the sensor at the hot end and the second at the cold end
		heatPoints is a list of two space indexes used to determine the heat flux for the hot end (first term) and
			cold end (second term)
		mode is the mode used for the power calculations (e.g. COP) performed at the end of the simulation
		"""



		#restingTimes definition
		if restingTimeHot=='default':
			restingTimeHot=0.

		if restingTimeCold=='default':
			restingTimeCold=0.


		# information for the log file
		print ''
		print ''
		print '######################################################'
		print ''
		print '------------------------------------------------------'
		print fileName
		print '------------------------------------------------------'
		print ''
		print 'heatconpy version:',version
		print 'Module: magcalsys_solidstate_1D'
		if note != None:
			print ''
			print 'Note:',note
		print ''
		print 'Mode:', mode
		print 'System:', leftReservoir_material+'/'+leftThermalSwitch_material+'/'+MCM_material+'/'+rightThermalSwitch_material+'/'+rightReservoir_material
		print 'Dimensions (m):', str(dx*leftReservoir_length)+'/'+str(dx*leftThermalSwitch_length)+'/'+str(dx*MCM_length)+'/'+str(dx*rightThermalSwitch_length)+'/'+str(dx*rightReservoir_length)
		print 'Number of points:', str(leftReservoir_length)+'/'+str(leftThermalSwitch_length)+'/'+str(MCM_length)+'/'+str(rightThermalSwitch_length)+'/'+str(rightReservoir_length)
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

		#used to calculate the simulation time at the end
		startTime=time.time()

		if heatPoints=='default':
			leftHeatSensor=leftReservoir_length-3
			rightHeatSensor=-rightReservoir_length+3
		else:
			leftHeatSensor=heatPoints[0]
			rightHeatSensor=heatPoints[1]		

		if startingField != 'magnetization':
			initialState=True
		else:
			initialState=False

		#initializes the object for the simulation
		a=heatcond.heatcond_activemat_1D(ambTemperature, dx=dx, dt=dt, fileName=fileName, 
			materials=[leftReservoir_material,leftThermalSwitch_material,MCM_material,rightThermalSwitch_material,rightReservoir_material],
			borders=[1,leftReservoir_length+1,leftReservoir_length+leftThermalSwitch_length+1,leftReservoir_length+leftThermalSwitch_length+MCM_length+1,leftReservoir_length+leftThermalSwitch_length+MCM_length+rightThermalSwitch_length+1,leftReservoir_length+leftThermalSwitch_length+MCM_length+rightThermalSwitch_length+rightReservoir_length+1], 
			materialsOrder=[0,1,2,3,4],boundaries=boundaries,heatPoints=[leftHeatSensor,rightHeatSensor],initialState=initialState)

		#defines some variable for the cycles
		if temperatureSensor=='default':
			rightTemperatureSensor=-(rightReservoir_length/2)
			leftTemperatureSensor=(leftReservoir_length/2)
		else:
			rightTemperatureSensor=temperatureSensor[1]
			leftTemperatureSensor=temperatureSensor[0]

		value1=ambTemperature
		value2=ambTemperature
		i=0
		period=(1./freq)
		stopCriteria2=0.
		maximumPower=0.
		maximumWorkingPower = 0.
		maximumCOP = 0.

		#cycle simulation
		while (abs((value1-value2)/value2)>stopCriteria or i<minCycleNumber or abs((value1-value2)/value2)>=stopCriteria2) and i<maxCycleNumber:

			stopCriteria2 = abs((value1-value2)/value2)
			heatLeft=a.heatLeft
			heatRight=a.heatRight


			if i !=0 or startingField=='magnetization':

			#magnetization and computation

				#mode 1
				if magnetizationMode == 'constant_right':
					for j in range(magnetizationSteps):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						a.compute(((1-restingTimeHot-restingTimeCold)*period/2.)/magnetizationSteps,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				#mode 2
				if magnetizationMode == 'accelerated_right':
					for j in range(magnetizationSteps):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						#delta_t = (1/(2*freq*np.sqrt(magnetizationSteps)))*(np.sqrt(j+1)-np.sqrt(j))
						delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(magnetizationSteps)))*(np.sqrt(j+1)-np.sqrt(j))
						a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				#mode 3
				if magnetizationMode == 'decelerated_right':
					for j in range(magnetizationSteps):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(magnetizationSteps)))*(np.sqrt(magnetizationSteps-j)-np.sqrt(magnetizationSteps-j-1))
						a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				#mode 4
				if magnetizationMode == 'constant_left':
					for j in range(magnetizationSteps-1,-1,-1):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						a.compute(((1-restingTimeHot-restingTimeCold)*period/2.)/magnetizationSteps,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				#mode 5
				if magnetizationMode == 'accelerated_left':
					for j in range(magnetizationSteps-1,-1,-1):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(magnetizationSteps)))*(np.sqrt(magnetizationSteps-j)-np.sqrt(magnetizationSteps-j-1))
						a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				#mode 6
				if magnetizationMode == 'decelerated_left':
					for j in range(magnetizationSteps-1,-1,-1):
						a.activate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/magnetizationSteps+1)
						delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(magnetizationSteps)))*(np.sqrt(j+1)-np.sqrt(j))
						a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

				if restingTimeHot != 0.:
					a.compute(restingTimeHot*period,int(1/(freq*dt*cyclePoints)),solver=solverMode)




			#magnetization and computation

			#mode 1
			if demagnetizationMode == 'constant_right':
				for j in range(demagnetizationSteps):
					a.deactivate(leftReservoir_length+ j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					a.compute(((1-restingTimeHot-restingTimeCold)*period/2.)/demagnetizationSteps,int(1/(freq*dt*cyclePoints)),solver=solverMode)

			#mode 2
			if demagnetizationMode == 'accelerated_right':
				for j in range(demagnetizationSteps):
					a.deactivate(leftReservoir_length+ j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(demagnetizationSteps)))*(np.sqrt(j+1)-np.sqrt(j))
					a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

			#mode 3
			if demagnetizationMode == 'decelerated_right':
				for j in range(demagnetizationSteps):
					a.deactivate(leftReservoir_length+ j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(demagnetizationSteps)))*(np.sqrt(demagnetizationSteps-j)-np.sqrt(demagnetizationSteps-j-1))
					a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

			#mode 4
			if demagnetizationMode == 'constant_left':
				for j in range(demagnetizationSteps-1,-1,-1):
					a.deactivate(leftReservoir_length+ j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					a.compute(((1-restingTimeHot-restingTimeCold)*period/2.)/demagnetizationSteps,int(1/(freq*dt*cyclePoints)),solver=solverMode)

			#mode 5
			if demagnetizationMode == 'accelerated_left':
				for j in range(demagnetizationSteps-1,-1,-1):
					a.dectivate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(demagnetizationSteps)))*(np.sqrt(demagnetizationSteps-j)-np.sqrt(demagnetizationSteps-j-1))
					a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)

			#mode 6
			if demagnetizationMode == 'decelerated_left':
				for j in range(demagnetizationSteps-1,-1,-1):
					a.deactivate(leftReservoir_length+j*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1,leftReservoir_length+(j+1)*(leftThermalSwitch_length+MCM_length+rightThermalSwitch_length)/demagnetizationSteps+1)
					delta_t = (1/(2*(1./(1.-restingTimeHot-restingTimeCold))*freq*np.sqrt(demagnetizationSteps)))*(np.sqrt(j+1)-np.sqrt(j))
					a.compute(delta_t,int(1/(freq*dt*cyclePoints)),solver=solverMode)


			if mode=='heat_pump':
				if -(a.heatLeft-heatLeft)>maximumPower:
					maximumPower=-(a.heatLeft-heatLeft)/period
					maximumWorkingPower = -((a.heatLeft-heatLeft)-(a.heatRight-heatRight))/period
					maximumCOP = -(a.heatLeft-heatLeft)/(-(a.heatLeft-heatLeft)+(a.heatRight-heatRight))

			if mode=='refrigerator':
				if -(a.heatRight-heatRight)>maximumPower:
					maximumPower=-(a.heatRight-heatRight)/period
					maximumWorkingPower = -((a.heatLeft-heatLeft)-(a.heatRight-heatRight))/period
					maximumCOP = -(a.heatRight-heatRight)/(-(a.heatLeft-heatLeft)+(a.heatRight-heatRight))

			if restingTimeCold != 0.:
				a.compute(restingTimeCold*period,int((((1-restingTimeHot-restingTimeCold)*period)/dt)/(cyclePoints/2)),solver=solverMode)

			#updates the error values and prints information for the log file
			value1=value2
			value2=a.temperature[-5][1]
			i=1+i

		#prints more information for the log file
		endTime =time.time()
		simulationTime = endTime-startTime
		hours = int(simulationTime / 3600)
		minutes = int((simulationTime - hours*3600)/ 60)
		seconds = int(simulationTime - hours*3600 - (minutes * 60))
		hours = '%02d' % hours
		minutes = '%02d' % minutes
		seconds = '%02d' % seconds

		print '------------------------------------------------------'
		print ''
		print 'Number of cycles:', i
		print 'Final cycle error:', abs((value1-value2)/value2)
		print 'Maximum thermal power (W):', maximumPower
		print 'Maximum working power (W)', maximumWorkingPower
		print 'Maximum COP:', maximumCOP
		print 'No load temperature span (K):', -a.temperature[rightTemperatureSensor][1]+a.temperature[leftTemperatureSensor][1]
		print 'Working power at no load temperature span (W):', -((a.heatLeft-heatLeft)-(a.heatRight-heatRight))/period
		if mode=='heat_pump':
			print 'COP at no load temperature span:', -(a.heatLeft-heatLeft)/(-(a.heatLeft-heatLeft)+(a.heatRight-heatRight))
		if mode=='refrigerator':
			print 'COP at no load temperature span:', -(a.heatRight-heatRight)/(-(a.heatLeft-heatLeft)+(a.heatRight-heatRight))
		print 'Final time (s):', a.timePassed
		print 'Simulation duration:', hours+':'+minutes+':'+seconds
		print ''
		print '------------------------------------------------------'
		print ''
		print ''
		print ''