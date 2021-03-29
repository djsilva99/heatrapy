"""dimension_1.

This submodule allows to create unidimensional thermal objects, establish
thermal contact between them, activate or deactivate the whole, or part, of the
materials, add heat sources, modify thermal properties, define boundary
conditions, and compute the respective heat transfer processes.

The submodule is based on three classes:

Object
--------
This class only creates a single unidimensional thermal object. It includes
several methods that manipulate the thermal object.

SingleObject
--------------
This class computes the heat transfer processes involved in only one
unidimensional thermal object. It uses the class Object changing the state,
modifying the thermal properties and activate or deactivate the material.

SystemObjects
---------------
This class creates a system of unidimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object changing the state, modifying the thermal properties
and activate or deactivate the material.

"""

from .objects import Object, SystemObjects, SingleObject

__all__ = ['Object', 'SystemObjects', 'SingleObject']
