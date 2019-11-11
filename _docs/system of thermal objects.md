---
layout: home
permalink: /docs/system_of_thermal_objects/
author_profile: true
sidebar:
  nav: "docs"
---

The `system_objects` class can be used to compute heat transfer processes between solids. It creates a system of thermal objects, establishes contact between them and computes the respective thermal processes. To create a system of thermal objects type:

```python
>>> foo = ht.system_objects(
...     number_objects=2, materials=('Cu', 'Cu'),
...	objects_length=(10, 10), amb_temperature=293,
...	dx=0.01, dt=0.1, file_name='data',
...	initial_state=False,
...	boundaries=((2, 0), (3, 0))
...	)
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

The `contact` input parameter is a tuple of length 3 (one element for thermal object A, one for thermal object B, and one for the heat transfer coefficient). Each thermal object element is a tuple of length 2 where the first element is the index of the thermal object and the second is the spatial point index.

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
>>> foo.compute(
...     time_interval, write_interval,
...	solver='implicit_k')
...	)
```

This method computes the system for `time_interval`, and writes into the `file_name` file every `write_interval` time steps. Four different `solvers` can be used: `'explicit_general'`, `'explicit_k(x)'`, `'implicit_general'`, and `'implicit_k(x)'`. The implicit solvers use the Crank-Nicholsen method. While the solvers ending with general is only suited to x-independent thermal conductivities, the solvers ending with k(x) take into account x-dependent thermal conductivities. In general, the solver `implicit_k(x)` works for all the systems but is computationally more heavy (more time consuming).

