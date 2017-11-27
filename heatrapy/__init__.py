# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Heatrapy.

This module contains classes to extract the physical properties from
the materials from the heatrapy database, to create and compute
1 dimensional models of heat conduction, and to use the 1 dimensional
model to build models that include the magnetocaloric technology.

The package relies on 3 classes:

#########
calmatpro
#########

Extracts and interpolate the physical properties density (rho), specific
heat (cp), adiabatic change of temperature (tad), and thermal conductivity
(k), giving the temperature.

#########
heatcond
#########

Creates 1 dimensional models and computes the system using one of 4 solvers:
implicit with k(x), implicit with constant k, explicit with k(x), and explicit
with constant k. The output data is saved in a file.

#########
magcalsys
#########

Creates 1 dimensional models of magnetocaloric systems using the heatcond
class. During the computation, it writes log files of the key parameters that
can be analyzed at the end.

"""

from heatcond import heatcond_activemat_1D
from magcalsys import magcalsys_solidstate_1D

__all__ = [heatcond_activemat_1D, magcalsys_solidstate_1D]
