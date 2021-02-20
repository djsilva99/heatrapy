"""Contains the explicit_k solver.

Used to compute thermal processes

"""

import copy


def explicit_k(obj):
    """explicit_k solver.

    Used to compute one time step of systems with x-dependent thermal
    conductivity.

    """
    x = copy.deepcopy(obj.temperature)

    # computes
    for i in range(1, obj.num_points - 1):
        eta = obj.dt / (2. * obj.rho[i] * obj.Cp[i] * obj.dx * obj.dx)
        beta = obj.dt / (obj.rho[i] * obj.Cp[i])

        t_new = ((1 + beta * obj.Q[i]) * obj.temperature[i][0] +
                 eta * ((obj.k[i + 1] + obj.k[i]) * obj.temperature[i + 1][0] -
                        (obj.k[i - 1] + obj.k[i + 1] + 2 * obj.k[i]) *
                        obj.temperature[i][0] + (obj.k[i - 1] + obj.k[i]) *
                        obj.temperature[i - 1][0]) +
                 beta * (obj.Q0[i] - obj.Q[i] * obj.amb_temperature))

        x[i][1] = t_new

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

    # updates temperature for next iteration
    for i in range(0, obj.num_points):
        x[i][0] = x[i][1]

    nx = []
    for i in range(len(x)):
        nx.append(x[i][0])

    # latent heat
    lheat = copy.copy(obj.lheat)
    for i in range(1, obj.num_points - 1):
        j = 0
        for lh in obj.latent_heat[i]:
            temper = obj.temperature[i][0]
            if nx[i] > lh[0] and temper <= lh[0] and lheat[i][j][1] != lh[1]:
                en = obj.Cp[i] * obj.rho[i] * (nx[i] - obj.temperature[i][0])
                if en + lheat[i][j][1] >= lh[1]:
                    lheat[i][j][1] = lh[1]
                    energy_temp = lheat[i][j][1] + en - lh[1]
                    nx[i] = obj.temperature[i][0] + \
                        energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    nx[i] = obj.temperature[i][0]
            if nx[i] < lh[0] and temper >= lh[0] and lheat[i][j][1] != 0:
                en = obj.Cp[i] * obj.rho[i] * (nx[i] - obj.temperature[i][0])
                if en + lheat[i][j][1] <= 0.:
                    lheat[i][j][1] = 0.
                    energy_temp = (en + lheat[i][j][1])
                    nx[i] = obj.temperature[i][0] + \
                        energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    nx[i] = obj.temperature[i][0]
            j += 1

    y = copy.deepcopy(obj.temperature)

    # updates the temperature list
    for i in range(obj.num_points):
        y[i][1] = nx[i]
        y[i][0] = nx[i]

    return y, lheat
