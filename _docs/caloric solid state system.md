---
layout: home
permalink: /docs/caloric_solid_state_system/
author_profile: true
sidebar:
  nav: "docs"
---

To compute a fully solid state caloric system the function `solid_active_regenerator` must be called. The active regenerative processes can be used with the several allowed modes for application and removal of fields. Cascades of materials can also be computed. To run one simulation type

```python
>>> ht.solid_active_regenerator(
...     file_name, amb_temperature=293,
...	left_thermalswitch_length=2,
...	right_thermalswitch_length=2, MCM_length=20,
...	right_reservoir_length=3, left_reservoir_length=3,
...	MCM_material=((0.002, 'Gd'),),
...	left_thermalswitch_material='idealTS_hot',
...	right_thermalswitch_material='idealTS_cold',
...	left_reservoir_material='Cu',
...	right_reservoir_material='Cu', freq=.1, dt=.01,
...	dx=0.002, stop_criteria=5e-3,
...	solver='implicit_k(x)', min_cycle_number=30,
...	max_cycle_number=31, field_removal_steps=3,
...	field_applied_steps=1,
...	field_removal_mode='accelerated_right',
...	field_applied_mode='accelerated_left',
...	cycle_points=25, boundaries=(293, 293), note=None,
...	temperature_sensor='default',
...	heat_points='default', mode='refrigerator',
...	version=None, resting_time_hot='default',
...	resting_time_cold='default',
...	starting_field='applied',
...	type_study='fixed_temperature_span',
...	h_left=50000., h_right=50000., mod_freq='default',
...     materials_path=False
... )
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
* `materials_path`: string indicating the path to the materials folder. If False, then the materials forder is the builtin heatrapy database folder.
