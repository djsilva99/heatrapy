# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the explicit_general solver.

Used to compute thermal processes

"""

import numpy as np
import copy


def explicit_general(obj):
    """explicit_general solver.

    Used to compute one time step of systems with fixed thermal contuctivity.

    """

    x = copy.copy(obj.temperature)

    # computes
    for i in range(1, obj.num_points - 1):

        alpha = obj.dt * \
            obj.k[i] / (obj.rho[i] * obj.Cp[i] *
                        obj.dx * obj.dx)
        beta = obj.dt / (obj.rho[i] * obj.Cp[i])

        Tnew = ((1 + beta * obj.Q[i]) * obj.temperature[i][0] +
                alpha * (obj.temperature[i - 1][0] - 2 *
                         obj.temperature[i][0] + obj.temperature[i + 1][0]) +
                beta * (obj.Q0[i] - obj.Q[i] * obj.amb_temperature))
        x[i][1] = Tnew

    # left boundary for next time step
    if obj.boundaries[0] == 0:
        x[0][1] = obj.temperature[1][0]
    else:
        x[0][1] = obj.boundaries[0]

    # right boundary for next time step
    if obj.boundaries[1] == 0:
        x[obj.num_points - 1][1] = obj.temperature[obj.num_points - 2][0]
    else:
        x[obj.num_points - 1][1] = obj.boundaries[1]

    return x
