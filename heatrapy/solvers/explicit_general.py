import numpy as np
import copy


def explicit_general(obj):

    x=copy.copy(obj.temperature)

    # computes
    for i in range(1, obj.numPoints - 1):

        alpha = obj.dt * \
            obj.k[i] / (obj.rho[i] * obj.Cp[i] *
                        obj.dx * obj.dx)
        beta = obj.dt / (obj.rho[i] * obj.Cp[i])

        Tnew = (1 + beta * obj.Q[i]) * obj.temperature[i][0] \
            + alpha * (obj.temperature[i - 1][0] -
                    2 * obj.temperature[i][0] +
                    obj.temperature[i + 1][0]) + \
            beta * (obj.Q0[i] - obj.Q[i] * obj.ambTemperature)
        x[i][1] = Tnew

    # left boundary for next time step
    if obj.boundaries[0] == 0:
        x[0][1] = obj.temperature[1][0]
    else:
        x[0][1] = obj.boundaries[0]

    # right boundary for next time step
    if obj.boundaries[1] == 0:
        x[obj.numPoints - 1][1] =  obj.temperature[obj.numPoints - 2][0]
    else:
        x[obj.numPoints - 1][1] = obj.boundaries[1]

    return x
