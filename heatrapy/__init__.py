"""Heatrapy.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1D and 2D. It includes two system
models for the computation of caloric systems in 1D.

The package is based on three classes that create general models:

######
Object1D
######

This class only creates a single unidimensional thermal object. It includes two
methods: material activation and material deactivation, of part of the object.

######
Object2D
######

This class only creates a single two-dimensional thermal object. It includes
two methods: material activation and material deactivation, of part of the
object.

##############
SystemObjects1D
##############

This class creates a system of unidimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object1D for the creation of each thermal object.

##############
SystemObjects2D
##############

This class creates a system of two-dimensional thermal objects that can be in
contact to each other and computes the respective heat transfer processes. It
uses the class Object2D for the creation of each thermal object.

#############
SingleObject1D
#############

This class computes the heat transfer processes involved in only one
unidimensional thermal object. It uses the class Object1D for activating and
deactivating the material.

#############
SingleObject2D
#############

This class computes the heat transfer processes involved in only one
two-dimensional thermal object. It uses the class Object1D for activating and
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

from .dimension_1.objects import Object as Object1D
from .dimension_1.objects import SystemObjects as SystemObjects1D
from .dimension_1.objects import SingleObject as SingleObject1D
from .dimension_2.objects import Object as Object2D
from .dimension_2.objects import SystemObjects as SystemObjects2D
from .dimension_2.objects import SingleObject as SingleObject2D

__all__ = [Object1D, SystemObjects1D, SingleObject1D, Object2D,
           SystemObjects2D, SingleObject2D]
