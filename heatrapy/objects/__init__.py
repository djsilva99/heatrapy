# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Objects.

This submodule includes the object class for creating thermal objects,
object system class for creating systems of thermal objects, and
a complex single object for using in more complex thermal computing.

"""

from object import object
from system import system_objects, single_object

__all__ = [object, system_objects, single_object]
