"""Contains the 1D implicit_k solver.

Used to compute unidimensional thermal processes

"""

import numpy as np

from ._latent_heat import apply_latent_heat


def implicit_k(obj):
    """implicit_k solver.

    Used to compute one time step of 1D systems with k-dependent thermal
    conductivities.

    """
    n = obj.num_points

    # initializes the matrixes for the equation systems
    a = np.zeros((n, n))
    b = np.zeros(n)

    # left boundary
    a[0][0] = 1
    if obj.boundaries[0] == 0:
        b[0] = obj.temperature[1][0]
    else:
        b[0] = obj.boundaries[0]

    # right boundary
    a[n - 1][n - 1] = 1
    if obj.boundaries[1] == 0:
        b[n - 1] = obj.temperature[n - 2][0]
    else:
        b[n - 1] = obj.boundaries[1]

    # creates the matrixes and solves the equation systems
    for i in range(1, n - 1):
        gamma = 4. * obj.rho[i] * obj.Cp[i] * obj.dx * obj.dx / obj.dt

        a[i][i - 1] = obj.k[i - 1] + obj.k[i]
        a[i][i] = -(gamma + obj.k[i + 1] + obj.k[i - 1] +
                    2. * obj.k[i] - 2 * obj.dt * obj.dt *
                    obj.Q[i])
        a[i][i + 1] = obj.k[i + 1] + obj.k[i]
        b[i] = -(obj.k[i + 1] + obj.k[i]) * \
            obj.temperature[i + 1][0] + \
            (-gamma + obj.k[i + 1] + obj.k[i - 1] + 2. *
                obj.k[i] - 2 * obj.dt * obj.dt * obj.Q[i]) * \
            obj.temperature[i][0] - \
            (obj.k[i - 1] + obj.k[i]) * \
            obj.temperature[i - 1][0] - \
            4. * obj.dx * obj.dx * \
            (obj.Q0[i] - obj.Q[i] * obj.amb_temperature)

    x = np.linalg.solve(a, b)

    # latent heat
    nx_list = x.tolist()
    lheat = apply_latent_heat(nx_list, obj)

    # pack into [current, next] pairs expected by the caller
    y = [[nx_list[i], nx_list[i]] for i in range(n)]

    return y, lheat
