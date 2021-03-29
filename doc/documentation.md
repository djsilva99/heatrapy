# Documentation

Heatrapy stands for HEAt TRAnsfer in PYthon. It simulates dynamic 1D and 2D heat transfer processes in solids using the finite difference method. Heatrapy includes both the modeling of caloric effects and the incorporation of phase transitions.

author: Daniel Silva (djsilva@gmx.com) <br> current version: v2.0.1

Installation
------------
To install heatrapy use the pip package manager:
```bash
$ pip install heatrapy
```
To import the heatrapy module type in the python shell:
```python
>>> import heatrapy as htp
```

Thermal object
--------------
Thermal objects are the building block of the heatrapy module. Each thermal object is defined by a range of discrete space, where each point incorporates thermal properties, power sources and temperature. Boundary conditions are also defined and can be changed at any time. Thermal objects can use several methods that change the sate of each point, materials properties, and power sources. The state of each point defines if a field is applied (if True) or not (if False). This is of paramount importance for caloric materials. Moreover, phase change transitions can be incorporated according to latent heat data. There are two different types of thermal object classes: SingleObject and SystemObjects. While the first only computes heat transfer processes of single thermal objects, the second can compute heat transfer processes of several thermal objects that can contact each other.

Database of materials
---------------------
A builtin database of materials is provided in the database folder. This database is very limited. The use of customized data can be done by defining the material data absolute path variable `file_name`. The builtin database includes several standard materials and can be consulted <a href='https://github.com/djsilva99/heatrapy/tree/master/heatrapy/database'>here</a>. Some folder names end with `_mag`. This means that the state of each point defines if a magnetic field is applied.

Visualization
-------------
Although the temperature data can be saved for each point and time in a csv file, the package also allows live visualization when using the SingleObject classes. By default, when initializing the related objects, the plotting of the temperature begins. The different live plot types are defined in the `draw` parameter. So far, one can plot the temperature, materials, heat power sources and state for SingleObject2D. Only the temperature can be plotted when using SingleObject1D.

Example
-------
As an example, the code below computes 100 s of a Cu 2D solid with a water circle and a aluminum square.
```python
>>> import heatrapy as htp
>>> example = htp.SingleObject2D(
...     293,
...     material='Cu',
...     boundaries=(300, 0, 0, 0),
...     size=(20, 20)
... )
>>> example.change_material(
...     material='Al',
...     shape='square',
...     initial_point=(5, 2),
...     length=(4, 4)
... )
>>> example.change_material(
...     material='water',
...     shape='circle',
...     initial_point=(5, 13),
...     length=4
... )
>>> example.compute(100, 10)
```
The output at the end of the computation will be the following:
<img src="https://github.com/danieljosesilva/heatrapy/blob/master/img/example.png">

Changelog
---------
The changelog can be consulted <a href='https://github.com/djsilva99/heatrapy/tree/master/CHANGELOG.md'>here</a>.

License
-------
Heatrapy is licensed under the terms of the <a href='https://github.com/djsilva99/heatrapy/tree/master/LICENSE'>MIT License (MIT)</a>.

Contributing
------------
Any type of contributions are welcome. A road map for future releases can be consulted in the <a href='https://github.com/djsilva99/heatrapy/wiki'>wiki pages</a>. Besides contributing by coding, feel free to report bugs and suggest new features. When contributing, please first discuss the change you wish to make via issue or email (djsilva@gmx.com), and then test the new code according to the <a href='https://github.com/djsilva99/heatrapy/tree/master/CONTRIBUTING.md'> contributing file</a>.

Citing Heatrapy
-----------------
Please acknowledge heatrapy if it contributes to a project that leads to a publication by citing <a href='https://github.com/djsilva99/heatrapy/wiki'>this paper</a>.
