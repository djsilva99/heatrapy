# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the explicit_k solver.

Used to compute thermal processes

"""

import numpy as np
import copy


def explicit_k(obj):
    """explicit_k solver.

    Used to compute one time step of systems with x-dependent thermal
    contuctivity.

    """

    x = copy.copy(obj.temperature)

    # computes
    for i in range(1, obj.num_points - 1):

        eta = obj.dt / (2. * obj.rho[i] * obj.Cp[i] * obj.dx * obj.dx)
        beta = obj.dt / (obj.rho[i] * obj.Cp[i])

        Tnew = ((1 + beta * obj.Q[i]) * obj.temperature[i][0] +
                eta * ((obj.k[i + 1] + obj.k[i]) * obj.temperature[i + 1][0] -
                       (obj.k[i - 1] + obj.k[i + 1] + 2 * obj.k[i]) *
                       obj.temperature[i][0] + (obj.k[i - 1] + obj.k[i]) *
                       obj.temperature[i - 1][0]) +
                beta * (obj.Q0[i] - obj.Q[i] * obj.amb_temperature))

        x[i][1] = Tnew

    # left boundary for next time step
    if obj.boundaries[0] == 0:
        x[0][1] = obj.temperature[1][1]
    else:
        x[0][1] = obj.boundaries[0]

    # right boundary for next time step
    if obj.boundaries[1] == 0:
        x[obj.num_points - 1][1] = obj.temperature[obj.num_points - 2][1]
    else:
        x[obj.num_points - 1][1] = obj.boundaries[1]

    return x
