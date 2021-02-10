"""Objects.

This submodule includes the object class for creating thermal objects,
object system class for creating systems of thermal objects, and
a complex single object for using in more complex thermal computing.

"""

from .object import Object
from .system import SystemObjects, SingleObject

__all__ = [Object, SystemObjects, SingleObject]
