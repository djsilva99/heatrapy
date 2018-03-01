# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Heatrapy.

This module allows to create thermal objects, establish thermal contact between
them, activate or deactivate the whole, or part, of the materials, and compute
the respective heat transfer processes, in 1 dimension. It includes several
system models for the computation of several thermotechnologies, including
ferroic-based systems.

There are 3 classes that create general models:

######
object
######

This class only creates a single thermal object. It includes 2 methods:
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

At this moment there are 2 complex systems that can be computed:

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

from systems import solid_active_regenerator, fluid_active_regenerator
from objects import object, system_objects, single_object

__all__ = [object, system_objects, single_object, solid_active_regenerator,
           fluid_active_regenerator]
