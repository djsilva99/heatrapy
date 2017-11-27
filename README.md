# heatrapy

This package is a module for simulating heat transfer processes. At the moment, it is focused on heat conduction, and includes a subpackage for computing magnetocaloric systems. The heatcond_activemat_1D class can be used to simulate conduction of heat in 1 dimension.
To simulate magnetocaloric operating systems use the magcalsys class. For visualizing the produced data use the python library <a href='https://github.com/danieljosesilva/physplotlib'>physplotlib</a>.


author: Daniel Silva (djsilva99@gmail.com) <br> current version: v0.1.1 (alpha)

Validated features:
- [x] temperature computation of heatcond models
- [ ] computation of heatcond models with heat source
- [x] temperature computation of magcalsys models
- [ ] output powers and COPs for magcalsys models


## Table of contents

1. [Instalation](#instalation)
2. [Introduction](#introduction)
3. [heatcond module](#heatcond)
    1. [creation of the model](#heatcond-model)
    2. [state](#heatcond-state)
    3. [power source](#heatcond-power)
    4. [compute](#heatcond-compute)
    5. [example](#heatcond-example)
4. [magcalsys module](#magcalsys)


## 1. Instalation <a name="instalation"></a>

To install simply use the pip package manager:

```bash
$ pip install heatrapy
```

To import the heatrapy module type in the python shell:

```python
>>> import heatrapy as ht
```


## 2. Introduction <a name="introduction"></a>

This module contains classes to extract the physical properties from the materials from the heatrapy database, to create and compute 1 dimensional models of heat conduction, and to use the 1 dimensional model to build models that include the magnetocaloric technology. The package relies on 3 classes:

**calmatpro**

Extracts and interpolate the physical properties density (rho), specific heat (cp), adiabatic change of temperature (tad), and thermal conductivity (k), giving the temperature. The current database includes: Gd, Cu, idealTS_cold, idealTS_hot, and water.

**heatcond**

Creates 1-dimensional models and computes the system using one of 4 solvers: implicit with k(x), implicit with constant k, explicit with k(x), and explicit with constant k. The output data is saved in a file.

**magcalsys**

Creates 1-dimensional models of magnetocaloric systems using the heatcond class. During the computation, it writes log files of the key parameters that can be analyzed at the end.


## 3. heatcond <a name="heatcond"></a>

This class solves numerically the heat conduction equation in 1 dimension of active material(s). Three solvers can be used: explicit with x-independent k, explicit with x-dependent k, implicit with x-independent k, and implicit with x-dependent k. The class has 4 main methods: activate and deactivate, for the activation and deactivation of part of the solid, deactivate for the deactivation of part of the solid, changeHeatPower for changing the heat source, and compute for solving the equation for a given period of time. This class is suited for simulations involving caloric systems, e.g. magnetocaloric or electrocaloric systems.


#### i. creation of the model <a name="heatcond-model"></a>

To create a model called foo type:

```python
foo = ht.heatcond_activemat_1D(ambTemperature, materials=['Cu'], borders=[1, 11], materialsOrder=[0],
                               dx=0.01, dt=0.1, fileName='data.txt', boundaries=[0, 0], Q=[], Q0=[],
                               heatPoints=[1, -2], initialState=False)
```

The input variables are the following:

* `ambTemperature`: ambient temperature of the whole system
* `materials`: list of strings of all the used materials present in the materials database
* `borders`: list of the points where there is a change of material
* `materialsOrder`: list of the materials list indexes that defines the material seperated by borders
* `dx`: space step
* `dt`: times step
* `fileName`: file name where the temperature and heat flux are saved
* `boundaries`: list of two entries that define the boundary condition. If 0 the boundary condition is insulation
* `Q`: list of 3 entry lists that gives the fixed heat source coeficient. The first term is the initial space index where it is applied. The second is the final space index where it is applied. The third is the value of the coeficient.
* `Q0` is a list of 3 entry lists that gives the temperature dependent heat source coefficient. The first term is the initial space index where it is applied. The second is the final space index where it is applied. The third is the value of the coeficient.
* `heatPoints`: list of the space indexes that we want to extract the heat flux. Normally, the first term is the heat flux of the hot end and the second term is the heat flux of the cold end


#### ii. state <a name="heatcond-state"></a>

To activate a piece of material type:

```python
foo.activate(initialPoint, finalPoint)
```

This command will activate the material from the `initialPoint` to the `finalPoint`.

To deactivate a piece of material type:

```python
foo.deactivate(initialPoint, finalPoint)
```

This command will deactivate the material from the `initialPoint` to the `finalPoint`.


#### iii. power source <a name="heatcond-power"></a>

To change the initial power source type:

```python
foo.changeHeatPower(self, Q=[], Q0=[])
```

The `Q` and `Q0` variables are already described above.


#### iv. compute <a name="heatcond-compute"></a>

To compute the system type:

```python
foo.compute(timeInterval, writeInterval, solver='explicit_k(x)',
            modeTemp=False, numFlag=0.5, modeTempPoint=1)
```

This method computes the system for `timeInterval`, and writes into `fileName` every `writeInterval` time steps. Four different `solver` can be used: `'explicit_general'`, `'explicit_k(x)'`, `'implicit_general'`, and `'implicit_k(x)'`. If `modeTemp` is True, the `timeInterval` is no longer applied the computation stops when the temperature of `modeTempPoint` index changes `numFlag` relative to the initial value.


#### v. example <a name="heatcond-example"></a>

The following example computes a simple 1-dimensional model with 0.5 m of gadolinium (Gd). The system is initial at 293 K. One end of the system is at a fixed temperature of 300 K, while the other end is insulated. The used time step is 1 second, and the used space step is 0.05 m, so that the overal number of space points is 10. The system is initialy deactivated. To create the model we initialize the object `example`:

```python
example = heatcond.heatcond_activemat_1D(293, materials=['Gd'], borders=[1,11], 
                                         materialsOrder=[0], dx=0.05, dt=1.,
                                         fileName='example.txt', boundaries=[300,0],
                                         Q=[], Q0=[], initialState=False)
```

Then we compute the system for 30000 s, write the output values every 300 s, using the 'implicit_k(x)' solver:

```python
example.compute(30000,300,solver='implicit_k(x)')
```

Àfterwards we activate the whole system:

```python
example.activate(1,10)
```

and we compute the system for 30000 s one more time:

```python
example.compute(30000,300,solver='implicit_k(x)')
```

The output data is stored in file example.txt. To visualize the temperature as a function of time for the point index 3 we use the <a href='https://github.com/danieljosesilva/physplotlib'>physplotlib</a> library and write:

```python
import physplotlib as pp
example_visualization = pp.statplot()
example_visualization.loadFile('example.txt')
example_visualization.verticalPlot([0], [0], [[3]],y_title='temperature (K)')
```

This code will output the following time-dependent temperature plot:
<center>
<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/example.png" alt="Drawing" style="width: 80%;"/>
</center>
As expected the temperature will tend to the 300 K fixed temperature. The activation of the material at 30000 s makes the temperature to jump to ~ 302 K.

## 4. magcalsys <a name="magcalsys"></a>

The magcalsys module creates 1-dimensional models for magnetocaloric systems. At this moment it includes a class that computes 1-dimensional fully solid state magnetocaloric systems. The active magnetcic regenerative processes can be used with the several allowed demagnetization processes. Cascades of materials can also be computed.

The initialization and computation of the model is perfmed simultaneously by typing

```python
bar = ht.magcalsys_solidstate_1D(fileName, ambTemperature=293, leftThermalSwitch_length=10,
                                 rightThermalSwitch_length=10, MCM_length=50, rightReservoir_length=15,
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
```

The input variables are the following:

* `fileName`: name of the file where the temperature and heat flux is saved
* `ambTemperature`: ambient temperature of the whole system
* `leftThermalSwitch_length`: length of the left thermal switch
* `rightThermalSwitch_length`: length of the right thermal switch
* `leftThermalSwitch_material`: material of the left thermal switch
* `rightThermalSwitch_material`: material of the right thermal switch
* `MCM_length`: length of the magnetocaloric material
* `MCM_material`: magnetocaloric material
* `leftReservoir_length`: length of the left reservoir
* `rightReservoir_length`: length of the right reservoir
* `leftReservoir_material`: material of the left reservoir
* `rightReservoir_material`: material of the right reservoir
* `freq`: operating frequency
* `dt`: times step
* `dx`: space step
* `stopCriteria`: error threshold to stop the simulation
* `solverMode`: solver
* `minCycleNumber`: minimum number of cycles that has to be computed
* `maxCycleNumber`: maximum number of cycles that needs to be computed
* `demagnetizationSteps`: number of steps during the demagnetization
* `magnetizationSteps`: number of steps during the magnetization
* `demagnetizationMode`: mode of demagnetization. It can be constant_right, constant_left, accelerated_right, accelerated_left, decelerated_right, and decelerated_left
* `magnetizationMode` mode of magnetization: It can be constant_right, constant_left, accelerated_right, accelerated_left, decelerated_right, and decelerated_left
* `cyclePoints`: number of points recorded in each cycle
* `boundaries`: list with the boundary conditions used in heatcond
* `temperatureSensor`: list of two space indexes used to determine the temperature span at the end of the simulation. The first term is the sensor at the hot end and the second at the cold end
* `heatPoints`: list of two space indexes used to determine the heat flux for the hot end (first term) and cold end (second term)
* `mode`: mode used for the power calculations (e.g. COP) performed at the end of the simulation

