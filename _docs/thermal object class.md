---
layout: home
permalink: /docs/thermal_object_class/
author_profile: true
sidebar:
  nav: "docs"
---

1D thermal objects incorporates all thermal physical properties (temperature, specific heat, thermal conductivity and density) and boundary conditions, for a discretized length. Thermal objects can include heat sources and can be in contact with other thermal objects. There are three classes for thermal objects. The ```Object``` class is responsible for creating thermal objects, the ```SystemObjects``` class is responsible for creating systems of thermal objects and compute the respective system, and the ```SingleObject``` class is responsible for computing single objects made of several materials.

The `Object` class is the building block of the whole package. It creates thermal objects to be used in the more complex systems when the ```SystemObjects``` and ```SingleObject``` classes are called. It includes two methods to apply and remove fields. To create a thermal object type:

```python
>>> foo = ht.Object(
...     amb_temperature, materials=('Cu',), borders=(1, 11),
...	materials_order=(0,), dx=0.01, dt=0.1,
...	file_name=None, boundaries=(0, 0), Q=[], Q0=[],
... initial_state=False, heat_save=False,
... materials_path=False
... )
```

The input variables are the following:

* `amb_temperature`: ambient temperature of the whole system
* `materials`: tuple of strings of all the used materials present in the folder materials
* `borders`: tuple of the points where there is a change of material
* `materials_order`: tuple of the `materials` tuple index that defines the material properties given by `borders`
* `dx`: the space step
* `dt`: the times step
* `file_name`: file name where the temperature and heat flux are saved. If `None` no file is saved.
* `boundaries`: tuple of two entries that define the boundary condition for temperature. If 0 the boundary condition is insulation
* `Q`: list of fixed heat source coefficient.
* `Q0`: list of temperature dependent heat source coefficient.
* `initial_state`: initial state of the materials. `True` if the field is applied and `False` if the field is removed.
* `heat_save`: `True` if saving the heat at the two borders.
* `materials_path`: string indicating the path to the materials folder. If False, then the materials forder is the builtin heatrapy database folder.

In this case, the `foo` object will have the following main array attributes:
* `temperature`: temperature
* `k`: thermal conductivity
* `Cp`: specific heat
* `rho`: density
* `Q`: fixed heat source coefficient
* `Q0`: temperature dependent heat source coefficient

