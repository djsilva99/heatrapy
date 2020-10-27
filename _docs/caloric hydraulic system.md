---
layout: home
permalink: /docs/caloric_hydraulic_system/
author_profile: true
sidebar:
  nav: "docs"
---

To compute a hydraulic active caloric regenerative system the function `fluid_active_regenerator` must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type

```python
>>> ht.fluid_active_regenerator(
...     file_name, amb_temperature=298, fluid_length=160,
...	MCM_length=50, right_reservoir_length=20,
...	left_reservoir_length=20,
...	MCM_material=((0.000, 'Gd'),),
...	fluid_material='water',
...	left_reservoir_material='Cu',
...	right_reservoir_material='Cu', freq=0.17, dt=0.1,
...	dx=0.001, stop_criteria=1e-6,
...	solver='implicit_k(x)', min_cycle_number=1,
...	max_cycle_number=5000, cycle_points=25, note=None,
...	boundaries=((2, 298), (3, 0)), mode='heat_pump',
...	version=None, leftHEXpositions=10,
...	rightHEXpositions=10, starting_field='applied',
...	temperature_sensor=(3, -3), field_removal_steps=1,
...	field_applied_steps=1,
...	field_removal_mode='accelerated_left',
...	field_applied_mode='accelerated_right',
...	applied_static_field_time_ratio=(0., 0.),
...	removed_static_field_time_ratio=(0., 0.),
...	h_mcm_fluid=1, h_leftreservoir_fluid=1,
...	h_rightreservoir_fluid=1,
...	mcm_discontinuity='default', type_study='no_load',
...	stroke=.02, mod_freq='default', materials_path=False
... )
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
* `materials_path`: string indicating the path to the materials folder. If False, then the materials forder is the builtin heatrapy database folder.
