"""2D Objects.

This submodule includes the object class for creating two-dimensional thermal
objects, object system class for creating systems of two-dimensional thermal
objects, and a complex single object for using in more complex thermal
computing.

"""

from .object import Object
from .system import SystemObjects
from .single import SingleObject

__all__ = [Object, SystemObjects, SingleObject]
