"""dimension_1.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1D. It includes two system
models for the computation of caloric systems in 1D.

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


At this moment there are two type of caloric systems (functions) that can be
computed:

########################
fluid_active_regenerator
########################

This function creates and computes a unidimensional active regenerative system
used for refrigeration and heat pumps. The heat exchanger is a fluid. It can be
used to compute caloric systems, e.g. magnetocaloric, electrocaloric,
elastocaloric, and barocaloric.

########################
solid_active_regenerator
########################

This function creates and computes a unidimensional active regenerative system
used for refrigeration and heat pumps. The heat exchanger is the solid material
itself. It can be used to compute caloric systems, e.g. magnetocaloric,
electrocaloric, elastocaloric, and barocaloric.

"""

from .systems import solid_active_regenerator, fluid_active_regenerator
from .objects import Object, SystemObjects, SingleObject

__all__ = [Object, SystemObjects, SingleObject, solid_active_regenerator,
           fluid_active_regenerator]
