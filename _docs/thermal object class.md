---
layout: home
permalink: /docs/thermal_object_class/
author_profile: true
sidebar:
  nav: "docs"
---

1D thermal objects incorporates all thermal physical properties (temperature, specific heat, thermal conductivity and density) and boundary conditions, for a discretized length. Thermal objects can include heat sources and can be in contact with other thermal objects. There are three classes for thermal objects. The ```object``` class is responsible for creating thermal objects, the ```system_objects``` class is responsible for creating systems of thermal objects and compute the respective system, and the ```single_object``` class is responsible for computing single objects made of several materials.

The `object` class is the building block of the whole package. It creates thermal objects to be used in the more complex systems when the ```system_objects``` and ```single_object``` classes are called. It includes two methods to apply and remove fields. To create a thermal object type:

```python
>>> foo = ht.object(
...        amb_temperature, materials=('Cu',), borders=(1, 11),
...		materials_order=(0,), dx=0.01, dt=0.1,
...		file_name='data.txt', boundaries=(0, 0), Q=[],
...		Q0=[], initial_state=False, heat_save=False
...	)
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

