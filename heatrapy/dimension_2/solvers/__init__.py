"""Solvers.

This submodule contains the solvers for the computation of the two-dimensional
thermal processes.

"""

from .explicit_k import explicit_k
from .explicit_general import explicit_general

__all__ = ['explicit_k', 'explicit_general']
