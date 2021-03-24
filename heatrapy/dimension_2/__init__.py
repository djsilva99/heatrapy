"""dimension_2.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 2D. It includes two system models
for the computation of caloric systems in 1D.

The package is based on three classes that create general models:

######
Object2D
######

This class only creates a single two-dimensional thermal object. It includes
two methods: material activation and material deactivation, of part of the
object.

##############
SystemObjects2D
##############

This class creates a system of two-dimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object2D for the creation of each thermal object.

#############
SingleObject2D
#############

This class computes the heat transfer processes involved in only one
two-dimensional thermal object. It uses the class Object1D for activating and
deactivating the material.

"""


from .objects import Object, SystemObjects, SingleObject

__all__ = [Object, SystemObjects, SingleObject]
