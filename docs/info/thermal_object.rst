Thermal object
==============

Thermal objects are the building block of the heatrapy module. Each
thermal object is defined by a range of discrete space, where each
point incorporates thermal properties, power sources and
temperature. Boundary conditions are also defined and can be changed
at any time. Thermal objects can use several methods that change the
sate of each point, materials properties, and power sources. The state
of each point defines if a field is applied (if True) or not (if
False). This is of paramount importance for caloric
materials. Moreover, phase change transitions can be incorporated
according to latent heat data. There are two different types of
thermal object classes: SingleObject and SystemObjects. While the
first only computes heat transfer processes of single thermal objects,
the second can compute heat transfer processes of several thermal
objects that can contact each other.
