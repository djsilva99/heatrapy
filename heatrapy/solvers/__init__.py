# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Solvers.

This submodule contains the solvers for the computation of thermal processes.

"""

from implicit_k import implicit_k
from implicit_general import implicit_general
from explicit_k import explicit_k
from explicit_general import explicit_general

__all__ = [implicit_k, implicit_general, explicit_k, explicit_general]
