import numpy as np
import copy


def implicit_k(obj):

    # initializes the matrixes for the equation systems
    a = np.zeros((obj.numPoints, obj.numPoints))
    b = np.zeros(obj.numPoints)

    # left boundary
    a[0][0] = 1
    if obj.boundaries[0] == 0:
        b[0] = obj.temperature[1][0]
    else:
        b[0] = obj.boundaries[0]

    # right boundary
    a[obj.numPoints - 1][obj.numPoints - 1] = 1
    if obj.boundaries[1] == 0:
        value = obj.temperature[obj.numPoints - 2][0]
        b[obj.numPoints - 1] = value
    else:
        b[obj.numPoints - 1] = obj.boundaries[1]

    # creates the matrixes and solves the equation systems
    for i in range(1, obj.numPoints - 1):
        gamma = 4. * obj.rho[i] * obj.Cp[i] * \
            obj.dx * obj.dx / obj.dt

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
            (obj.Q0[i] - obj.Q[i] * obj.ambTemperature)

    x = np.linalg.solve(a, b)

    y=copy.copy(obj.temperature)

    for i in range(obj.numPoints):
        y[i][1] = x[i]
        y[i][0] = x[i]

    return y
