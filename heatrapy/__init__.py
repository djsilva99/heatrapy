"""Heatrapy.

author: Daniel Silva (djsilva@gmx.com)
current version: v2.0.3
url: https://github.com/djsilva99/heatrapy

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1D and 2D.

The module is based on four main classes:

SingleObject1D
--------------
This class computes the heat transfer processes involved in only one
unidimensional thermal object. It uses the class Object1D changing the state,
modifying the thermal properties and activate or deactivate the material.

SingleObject2D
--------------
This class computes the heat transfer processes involved in only one
two-dimensional thermal object. It uses the class Object2D changing the state,
modifying the thermal properties and activate or deactivate the material.

SystemObjects1D
---------------
This class creates a system of unidimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object1D changing the state, modifying the thermal properties
and activate or deactivate the material.

SystemObjects2D
---------------
This class creates a system of two-dimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object2D changing the state, modifying the thermal properties
and activate or deactivate the material.

"""

from .dimension_1.objects import SystemObjects as SystemObjects1D
from .dimension_1.objects import SingleObject as SingleObject1D
from .dimension_2.objects import SystemObjects as SystemObjects2D
from .dimension_2.objects import SingleObject as SingleObject2D

__all__ = ['SystemObjects1D', 'SingleObject1D', 'SystemObjects2D',
           'SingleObject2D']
