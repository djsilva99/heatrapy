"""Contains the implicit_k solver.

Used to compute thermal processes

"""

import numpy as np
import copy


def implicit_k(obj):
    """implicit_k solver.

    Used to compute one time step of systems with k-dependent thermal
    contuctivity.

    """
    # initializes the matrixes for the equation systems
    a = np.zeros((obj.num_points, obj.num_points))
    b = np.zeros(obj.num_points)

    # left boundary
    a[0][0] = 1
    if obj.boundaries[0] == 0:
        b[0] = obj.temperature[1][0]
    else:
        b[0] = obj.boundaries[0]

    # right boundary
    a[obj.num_points - 1][obj.num_points - 1] = 1
    if obj.boundaries[1] == 0:
        value = obj.temperature[obj.num_points - 2][0]
        b[obj.num_points - 1] = value
    else:
        b[obj.num_points - 1] = obj.boundaries[1]

    # creates the matrixes and solves the equation systems
    for i in range(1, obj.num_points - 1):
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
    lheat = copy.copy(obj.lheat)
    for i in range(1, obj.num_points - 1):
        j = 0
        for lh in obj.latent_heat[i]:
            temper = obj.temperature[i][0]
            if x[i] > lh[0] and temper <= lh[0] and lheat[i][j][1] != lh[1]:
                en = obj.Cp[i] * obj.rho[i] * (x[i] - obj.temperature[i][0])
                if en + lheat[i][j][1] >= lh[1]:
                    lheat[i][j][1] = lh[1]
                    energy_temp = lheat[i][j][1] + en - lh[1]
                    x[i] = obj.temperature[i][0] + \
                        energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    x[i] = obj.temperature[i][0]
            if x[i] < lh[0] and temper >= lh[0] and lheat[i][j][1] != 0:
                en = obj.Cp[i] * obj.rho[i] * (x[i] - obj.temperature[i][0])
                if en + lheat[i][j][1] <= 0.:
                    lheat[i][j][1] = 0.
                    energy_temp = (en + lheat[i][j][1])
                    x[i] = obj.temperature[i][0] + \
                        energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    x[i] = obj.temperature[i][0]
            j += 1

    y = copy.deepcopy(obj.temperature)

    for i in range(obj.num_points):
        y[i][1] = x[i]
        y[i][0] = x[i]

    return y, lheat
