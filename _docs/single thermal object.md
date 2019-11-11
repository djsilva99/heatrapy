---
layout: home
permalink: /docs/single_thermal_object/
author_profile: true
sidebar:
  nav: "docs"
---

The `single_object` class solves numerically the heat conduction equation for a single thermal object, in which domains can have different materials. This class inherits the methods of the `object` class. To create a single thermal object type

```python
>>> foo = ht.single_object(
...     amb_temperature, materials=('Cu',),
...	borders=(1, 11), materials_order=(0,),
...	dx=0.01, dt=0.1, file_name='data.txt',
...	boundaries=(0, 0), Q=[], Q0=[],
...	heat_points=(1, -2),
...	initial_state=False, h_left=50000.,
...	h_right=50000.
... )
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
>>> foo.compute(
...     time_interval, write_interval,
...	solver='explicit_k(x)', modeTemp=False,
...	numFlag=0.5, modeTempPoint=1
... )
```

This method computes the system for `time_interval`, and writes into the `file_name` file every `write_interval` time steps. The four different solvers pointed out for the `system_objects` class can be used. `heat_points` is a list that defines the points where the heat flux are calculated if `modeTemp` is `True` the compute method stops when the point `modeTempPoint` changes `numFlag` relative to the initial value

