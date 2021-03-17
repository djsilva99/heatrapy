"""Heatrapy.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1D. It includes two system models
for the computation of caloric systems.

The package is based on three classes that create general models:

######
object
######

This class only creates a single thermal object. It includes two methods:
material activation and material deactivation, of part of the object.

##############
system_objects
##############

This class creates a system of objects that can be in contact to each other
and computes the respective heat transfer processes. It uses the class object
for the creation of each thermal object.

#############
single_object
#############

This class computes the heat transfer processes involved in only one thermal
object. It uses the class object for activating and deactivating the material.

At this moment there are two type of caloric systems (functions) that can be
computed:

########################
fluid_active_regenerator
########################

This function creates and computes an active regenerative system used for
refrigeration and heat pumps. The heat exchanger is a fluid. It can be used
to compute caloric systems, e.g. magnetocaloric, electrocaloric, elastocaloric,
and barocaloric.

########################
solid_active_regenerator
########################

This function creates and computes an active regenerative system used for
refrigeration and heat pumps. The heat exchanger is the solid material itself.
It can be used to compute caloric systems, e.g. magnetocaloric, electrocaloric,
elastocaloric, and barocaloric.

"""

from .dimension_1.systems import solid_active_regenerator as solid_active_regenerator_1D
from .dimension_1.systems import fluid_active_regenerator as fluid_active_regenerator_1D
from .dimension_1.objects import Object as Object1D
from .dimension_1.objects import SystemObjects as SystemObjects1D
from .dimension_1.objects import SingleObject as SingleObject1D
from .dimension_2.objects import Object as Object2D
from .dimension_2.objects import SystemObjects as SystemObjects2D
from .dimension_2.objects import SingleObject as SingleObject2D

#from .dimension_1.systems import solid_active_regenerator, fluid_active_regenerator
#from .dimension_1.objects import Object, SystemObjects, SingleObject

#__all__ = [Object, SystemObjects, SingleObject, solid_active_regenerator,
#           fluid_active_regenerator]

__all__ = [Object1D, SystemObjects1D, SingleObject1D, solid_active_regenerator_1D,
           fluid_active_regenerator_1D, Object2D]
