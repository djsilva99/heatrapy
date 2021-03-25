"""dimension_1.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1D.

The package is based on three classes that create general models:

######
Object
######

This class only creates a single unidimensional thermal object. It includes two
methods: material activation and material deactivation, of part of the object.

##############
SystemObjects
##############

This class creates a system of unidimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object1D for the creation of each thermal object.

#############
SingleObject
#############

This class computes the heat transfer processes involved in only one
unidimensional thermal object. It uses the class Object1D for activating and
deactivating the material.

"""

from .objects import Object, SystemObjects, SingleObject

__all__ = [Object, SystemObjects, SingleObject]
