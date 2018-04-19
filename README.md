# heatrapy

<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/heatrapy.png" alt="Drawing" height="30"/> heatrapy v0.2.10

This package is a module for simulating dynamic heat transfer processes involving caloric effects in 1.5D systems by using the finite difference method. It is focused on heat conduction, and includes two subpackages for computing caloric systems. For visualizing the output data use the python library <a href='https://github.com/danieljosesilva/physplotlib'>physplotlib</a>.

author: Daniel Silva (djsilva99@gmail.com) <br> current version: v0.2.10

![resSwitch-screenshot](https://github.com/danieljosesilva/heatrapy/blob/master/img/example.gif)

## Table of contents

1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Thermal objects](#thermal_objects)
    1. [object class](#object-class)
    2. [material state](#material-state)
    3. [system_objects class](#system-objects-class)
    4. [single_object class](#single-object-class)
    5. [example](#example)
4. [Caloric systems](#caloric)
    1. [solid state system](#solid-state)
    2. [hydraulic active regenerative system](#hydraulic)


## 1. Installation <a name="installation"></a>

To install heatrapy use the pip package manager:

```bash
$ pip install heatrapy
```

To import the heatrapy module type in the python shell:

```python
>>> import heatrapy as ht
```


## 2. Introduction <a name="introduction"></a>

This module allows to create thermal objects, establish thermal contact between them, activate or deactivate part of the materials, and compute the respective heat transfer processes, in 1D. The computation uses the finite difference method. It includes two system models for the computation of caloric systems.

The package is based on three classes that create generic models:

**object**

This class only creates a single thermal object. It includes two methods: material activation and material deactivation, of part of the object.

**system_objects**

This class creates a system of objects that can be in contact to each other and computes the respective heat transfer processes. It uses the class object for the creation of each thermal object.

**single_object**

This class computes the heat transfer processes involved in only one thermal object. It uses the class object for activating and deactivating the material.

At this moment there are two type of caloric systems (functions) that can be computed:

**solid_active_regenerator**

This function creates and computes an active regenerative system used for
refrigeration and heat pumps. The heat exchanger is the **solid material** itself. It can be used to compute caloric systems, e.g. magnetocaloric, electrocaloric, elastocaloric, and barocaloric.

**fluid_active_regenerator**

This function creates and computes an active regenerative system used for
magnetocaloric refrigeration and heat pumps. The heat exchanger is a **fluid**. It can be used to compute caloric systems, e.g. magnetocaloric, electrocaloric, elastocaloric, and barocaloric.


## 3. Thermal objects <a name="thermal_objects"></a>

1D thermal objects incorporates all thermal physical properties (temperature, specific heat, thermal conductivity and density) and boundary conditions, for a discretized length. Thermal objects can include heat sources and can be in contact with other thermal objects. There are three classes for thermal objects. The ```object``` class is responsible for creating thermal objects, the ```system_objects``` class is responsible for creating systems of thermal objects and compute the respective system, and the ```single_object``` class is responsible for computing single objects made of several materials.


### 1. Object class <a name="object-class"></a>

The `object` class is the building block of the whole package. It creates thermal objects to be used in the more complex systems when the ```system_objects``` and ```single_object``` classes are called. It includes two methods to apply and remove fields. To create a thermal object type:

```python
>>> foo = ht.object(amb_temperature, materials=('Cu',), borders=(1, 11),
...                 materials_order=(0,), dx=0.01, dt=0.1, file_name='data.txt',
...                 boundaries=(0, 0), Q=[], Q0=[], initial_state=False,
...                 heat_save=False)
```

The input variables are the following:

* `amb_temperature`: ambient temperature of the whole system
* `materials`: tuple of strings of all the used materials present in the folder materials
* `borders`: tuple of the points where there is a change of material
* `materials_order`: tuple of the `materials` tuple index that defines the material properties given by `borders`
* `dx`: the space step
* `dt`: the times step
* `file_name`: file name where the temperature and heat flux are saved
* `boundaries`: tuple of two entries that define the boundary condition for temperature. If 0 the boundary condition is insulation
* `Q`: list of fixed heat source coefficient.
* `Q0`: list of temperature dependent heat source coefficient.
* `initial_state`: initial state of the materials. `True` if the field is applied and `False` if the field is removed.
* `heat_save`: `True` if saving the heat at the two borders.

In this case, the `foo` object will have the following main array attributes:
* `temperature`: temperature
* `k`: thermal conductivity
* `Cp`: specific heat
* `rho`: density
* `Q`: fixed heat source coefficient
* `Q0`: temperature dependent heat source coefficient



### 2. material state <a name="material-state"></a>

Following the last subsection, to activate a piece of `foo` type:

```python
foo.activate(initial_point, final_point)
```

This command will activate the material from the `initial_point` to the `final_point`.

To deactivate a piece of material type:

```python
foo.deactivate(initial_point, final_point)
```

This command will deactivate the material from the `initial_point` to the `final_point`.


### 3. system_objects class <a name="system-objects"></a>

The `system_objects` class can be used to compute heat transfer processes between solids. It creates a system of thermal objects, establishes contact between them and computes the respective thermal processes. To create a system of thermal objects type:

```python
>>> foo = ht.system_objects(number_objects=2, materials=('Cu', 'Cu'),
...                         objects_length=(10, 10), amb_temperature=293, dx=0.01,
...                         dt=0.1, file_name='data', initial_state=False,
...                         boundaries=((2, 0), (3, 0)))
```

The input variables are the following:

* `amb_temperature`: ambient temperature of the whole system
* `materials`: tuple of strings of all the used materials present in the folder materials
* `number_objects`: integer for the number of thermal objects
* `objects_length`: tuple of the object lengths (spacial steps)
* `dx`: the space step
* `dt`: the times step
* `file_name`: file name where the temperature and heat flux are saved
* `boundaries`: tuple of two entries that define the boundary condition for temperature. The first corresponds to the thermal object while the second defines the temperature. If 0 the boundary condition is insulation
* `initial_state`: initial state of the materials. True if applied field and False is removed field.

To add a contact between thermal objects type

```python
>>> foo.contactAdd(contact)
```

The `contact` input parameter is a tuple of length 2 (one for thermal object A and one for thermal object B) in which each element is a tuple of length 2 as well where the first element is the index of the thermal object and the second is the index spatial point.

To remove a specific `contact` type

```python
>>> foo.contactRemove(contact)
```

To filer all contacts of a material type

```python
>>> foo.contactFilter(object)
```

where `object` is the thermal object index.

Finally, to compute the overall system type

```python
>>> foo.compute(timeInterval, writeInterval, solver='implicit_k'))
```

This method computes the system for `timeInterval`, and writes into the `file_name` file every `writeInterval` time steps. Four different `solvers` can be used: `'explicit_general'`, `'explicit_k(x)'`, `'implicit_general'`, and `'implicit_k(x)'`. The implicit solvers use the Crank-Nicholsen method. While the solvers ending with general is only suited to x-independent thermal conductivities, the solvers ending with k(x) take into account x-dependent thermal conductivities. In general, the solver `implicit_k(x)` works for all the systems but is computationally more heavy (more time consuming)


### 4. single_object class <a name="single-object"></a>

The `single_object` class solves numerically the heat conduction equation for a single thermal object, in which domains can have different materials. This class inherits the methods of the `object` class. To create a single thermal object type

```python
>>> foo = ht.single_object(amb_temperature, materials=('Cu',), borders=(1, 11),
...                        materials_order=(0,), dx=0.01, dt=0.1, file_name='data.txt',
...                        boundaries=(0, 0), Q=[], Q0=[], heat_points=(1, -2),
...                        initial_state=False, h_left=50000., h_right=50000.)
```

The input variables are the following:

* `amb_temperature`: ambient temperature of the whole system
* `materials`: tuple of strings of all the used materials present in the folder materials
* `borders`: tuple of the points where there is a change of material
* `materials_order`: tuple of the `materials` list indexes that defines the material properties given by `borders`
* `dx`: the space step
* `dt`: the times step
* `file_name`: file name where the temperature and heat flux are saved
* `boundaries`: tuple of two entries that define the boundary condition for temperature. If 0 the boundary condition is insulation
* `Q`: list of fixed heat source coefficient.
* `Q0`: list of temperature dependent heat source coefficient.
* `heat_points`: tuple of the space indexes where we want to extract the heat flux. Normally, the first term is the heat flux of the hot end and the second term is the heat flux of the cold end
* `initial_state`: initial state of the materials. `True` if the field is applied and `False` if the field is removed.
* `h_left`: left heat transfer coefficient
* `h_right`: right heat transfer coefficient

To compute the `single_object` type

```python
>>> foo.compute(timeInterval, writeInterval, solver='explicit_k(x)',
...             modeTemp=False, numFlag=0.5, modeTempPoint=1):
```

This method computes the system for `timeInterval`, and writes into the `file_name` file every `writeInterval` time steps. The four different solvers pointed out for the `system_objects` class can be used. `heat_points` is a list that defines the points where the heat flux are calculated if `modeTemp` is `True` the compute method stops when the point `modeTempPoint` changes `numFlag` relative to the initial value


### 5. example <a name="example"></a>

The following example computes a simple 1-dimensional model with 0.5 m of gadolinium (Gd). The system is initial at 293 K. One end of the system is at a fixed temperature of 300 K, while the other end is insulated. The used time step is 1 second, and the used space step is 0.05 m, so that the overall number of space points is 10. The system is initially deactivated. To create the model we initialize the object `example`:

```python
>>> example = ht.single_object(293, materials=['Gd'], borders=[1,11], materialsOrder=[0],
...                            dx=0.05, dt=1., fileName='example.txt', boundaries=[300,0],
...                            Q=[], Q0=[], initialState=False)
```

Then we compute the system for 30000 s, write the output values every 300 s, using the 'implicit_k(x)' solver:

```python
>>> example.compute(30000, 300, solver='implicit_k(x)')
```

Afterwards we activate the whole system:

```python
example.activate(1, 10)
```

and we compute the system for 30000 s one more time:

```python
>>> example.compute(30000, 300, solver='implicit_k(x)')
```

The output data is stored in file example.txt. To visualize the temperature as a function of time for the point index 3 we use the <a href='https://github.com/danieljosesilva/physplotlib'>physplotlib</a> library and write:

```python
>>> import physplotlib as pp
>>> example_visualization = pp.statplot()
>>> example_visualization.loadFile('example.txt')
>>> example_visualization.verticalPlot([0], [0], [[3]], y_title='temperature (K)')
```

This code will output the following time-dependent temperature plot:
<center>
<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/example.png" alt="Drawing" style="width: 80%;"/>
</center>
As expected the temperature will tend to the 300 K fixed temperature. The activation of the material at 30000 s makes the temperature to jump to ~ 302 K.


## 4. Caloric systems <a name="caloric"></a>

The purpose of the heatrapy package is to provide a framework in computing heat transfer processes in solids involving caloric effects, so that different systems can be computed. The developed model systems are included in the systems folder. At the momento it includes two different systems for caloric devices: a fully solid state device and an hydraulic active regenerative device. The computation is performed by executing the respective function, which is described for both systems below.


## 1. solid state system <a name="solid-state"></a>

To compute a fully solid state caloric system the function `solid_active_regenerator` must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type

```python
>>> ht.solid_active_regenerator(file_name, amb_temperature=293,
...                             left_thermalswitch_length=2,
...                             right_thermalswitch_length=2, MCM_length=20,
...                             right_reservoir_length=3, left_reservoir_length=3,
...                             MCM_material=((0.002, 'Gd'),),
...                             left_thermalswitch_material='idealTS_hot',
...                             right_thermalswitch_material='idealTS_cold',
...                             left_reservoir_material='Cu',
...                             right_reservoir_material='Cu',
...                             freq=.1, dt=.01, dx=0.002, stop_criteria=5e-3,
...                             solver='implicit_k(x)', min_cycle_number=30,
...                             max_cycle_number=31, field_removal_steps=3,
...                             field_applied_steps=1,
...                             field_removal_mode='accelerated_right',
...                             field_applied_mode='accelerated_left',
...                             cycle_points=25, boundaries=(293, 293), note=None,
...                             temperature_sensor='default',
...                             heat_points='default', mode='refrigerator',
...                             version=None, resting_time_hot='default',
...                             resting_time_cold='default',
...                             starting_field='applied',
...                             type_study='fixed_temperature_span',
...                             h_left=50000., h_right=50000.,
...                             mod_freq='default')
```

The input variables are the following:

* `file_name`: file name where the temperature and heat flux are saved
* `amb_temperature`: ambient temperature of the whole system
* `left_thermalswitch_length`: length of the left thermal switch
* `right_thermalswitch_length`: length of the right thermal switch
* `left_thermalswitch_material`: string for the material of the left thermal switch
* `right_thermalswitch_material`: string for the material of the right thermal switch
* `MCM_length`: length of the caloric material
* `MCM_material`: string for the material of the caloric material
* `left_reservoir_length`: length of the left reservoir
* `right_reservoir_length`: length of the right reservoir
* `left_reservoir_material`: string for the material of the left reservoir
* `right_reservoir_material`: string for the material of the right reservoir
* `freq`: operating frequency
* `dt`: times step
* `dx`: space step
* `stop_criteria`: error threshold to stop the simulation
* `solver`: solver
* `min_cycle_number`: minimum number of cycles that has to be computed
* `max_cycle_number`: maximum number of cycles that has to be computed
* `field_removal_steps`: number of steps during the field removal
* `field_applied_steps`: number of steps during the application of field
* `field_removal_mode`: mode of field removal modes can be `constant_right`, `constant_left`, `accelerated_right`, `accelerated_left`, `decelerated_right`, and `decelerated_left`
* `field_applied_mode`: the mode of the application of field modes can be `constant_right`, `constant_left`, `accelerated_right`, `accelerated_left`, `decelerated_right`, and `decelerated_left`
* `cycle_points`: number of points recorded for each position for each cycle
* `boundaries`: tuple with the boundary conditions
* `temperature_sensor`: tuple of two space indexes used to determine the temperature span at the end of the simulation. The first term is the sensor at the hot end and the second at the cold end
* `heat_points`: tuple of two space indexes used to determine the heat flux for the hot end (first term) and cold end (second term)
* `mode`: mode used for the power calculations (e.g. COP) performed at the end of the simulation. It can be `'refrigerator'` or `'heat_pump'`
* `version`: heatrapy version (default is None)
* `type_study`: 'no_load' or 'fixed_temperature_span'
* `h_left`: left heat transfer coefficient
* `h_right`: right heat transfer coefficient
* `mod_freq`: if not `'default'`, i.e. if tuple, allows to modulate the frequency according to a specific temperature. The first element is the file_name, and second the sensor point.


## 2. hydraulic active regenerative system <a name="hydraulic"></a>

To compute a hydraulic active caloric regenerative system the function `fluid_active_regenerator` must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type

```python
>>> ht.fluid_active_regenerator(file_name, amb_temperature=298, fluid_length=160,
...                             MCM_length=50, right_reservoir_length=20,
...                             left_reservoir_length=20,
...                             MCM_material=((0.000, 'Gd'),),
...                             fluid_material='water',
...                             left_reservoir_material='Cu',
...                             right_reservoir_material='Cu', freq=0.17, dt=0.1,
...                             dx=0.001, stop_criteria=1e-6,
...                             solver='implicit_k(x)', min_cycle_number=1,
...                             max_cycle_number=5000, cycle_points=25, note=None,
...                             boundaries=((2, 298), (3, 0)), mode='heat_pump',
...                             version=None, leftHEXpositions=10,
...                             rightHEXpositions=10, starting_field='applied',
...                             temperature_sensor=(3, -3), field_removal_steps=1,
...                             field_applied_steps=1,
...                             field_removal_mode='accelerated_left',
...                             field_applied_mode='accelerated_right',
...                             applied_static_field_time_ratio=(0., 0.),
...                             removed_static_field_time_ratio=(0., 0.),
...                             h_mcm_fluid=1,
...                             h_leftreservoir_fluid=1,
...                             h_rightreservoir_fluid=1,
...                             mcm_discontinuity='default',
...                             type_study='no_load', stroke=.02,
...                             mod_freq='default'):
```

The input variables are the following:

* `file_name`: file name where the temperature and heat flux are saved
* `amb_temperature`: ambient temperature of the whole system
* `fluid_length`: length of the fluid
* `fluid_material`: string for the material of the fluid
* `MCM_length`: length of the magnetocaloric material
* `MCM_material`: string for the material of the magnetocaloric material
* `left_reservoir_length`: length of the left reservoir
* `right_reservoir_length`: length of the right reservoir
* `left_reservoir_material`: string for the material of the left reservoir
* `right_reservoir_material`: string for the material of the right reservoir
* `freq`: operating frequency
* `dt`: times step
* `dx`: space step
* `stop_criteria`: error threshold to stop the simulation
* `solver`: solver
* `min_cycle_number`: minimum number of cycles that has to be computed
* `max_cycle_number`: maximum number of cycles that has to be computed
* `field_removal_steps`: number of steps during the field removal
* `field_applied_steps`: number of steps during the application of field
* `field_removal_mode`: mode of the field removal modes can be `constant_right`, `constant_left`, `accelerated_right`, `accelerated_left`, `decelerated_right`, and `decelerated_left`
* `field_applied_mode` is the mode of the application of field modes can be `constant_right`, `constant_left`, `accelerated_right`, `accelerated_left`, `decelerated_right`, and `decelerated_left`
* `cycle_points`: number of points recorded for each position for each cycle
* `boundaries`: tuple of two entries that define the boundary condition for temperature. The first corresponds to the thermal object while the second defines the temperature. If 0 the boundary condition is insulation
* `temperature_sensor`: tuple of two space indexes used to determine the temperature span at the end of the simulation. The first term is the sensor at the hot end and the second at the cold end
* `mode`: mode used for the power calculations (e.g. COP) performed at the end of the simulation. It can be `'refrigerator'` or `'heat_pump'`
* `version`: heatrapy version (default is None)
* `type_study`: `'no_load'` or `'fixed_temperature_span'`
* `h_mcm_fluid`: heat transfer coefficient for fluid - MCM
* `h_leftreservoir_fluid`: heat transfer coefficient for fluid - left reservoir
* `h_rigthreservoir_fluid`: heat transfer coefficient for fluid - right reservoir
* `mod_freq`: if not `'default'`, i.e. if a tuple, it allows to modulate the frequency according to a specific temperature. The first element of the tuple is the file_name, and second the sensor point.
* `velocity`: velocity of the fluid
* `leftHEXpositions`: distance in points from the left reservoir to the MCM
* `rightHEXpositions`: distance in points from the right reservoir to the MCM.
* `applied_static_field_time_ratio`: tuple with the resting time ratios before and after the application of the field.
* `removed_static_field_time_ratio`: tuple with the resting time ratios before and after the removal of the field.
* `mcm_discontinuity`: if not `'default'` the MCM is divided into n pieces. The input is a tuple where the first entry is the number of discontinuities, while the second is the thickness of each discontinuity in meters.
