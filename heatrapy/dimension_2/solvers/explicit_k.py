"""Contains the 2D explicit_k solver.

Used to compute two-dimensional thermal processes

"""

import copy


def explicit_k(obj):
    """explicit_k solver.

    Used to compute one time step of 2D systems with position dependent thermal
    conductivities.

    """
    x = copy.deepcopy(obj.temperature)
    temp_k = copy.deepcopy(obj.k)
    # bottom boundary
    x.append([])
    temp_k.append([])
    if obj.boundaries[1] == 0:
        for j in range(obj.size[1]):
            x[-1].append(obj.temperature[-1][j])
            temp_k[-1].append(obj.k[-1][j])
    else:
        for j in range(obj.size[1]):
            x[-1].append([obj.boundaries[1], obj.boundaries[1]])
            temp_k[-1].append(obj.k[-1][j])

    # top boundary
    x.insert(0, [])
    temp_k.insert(0, [])
    if obj.boundaries[0] == 0:
        for j in range(obj.size[1]):
            x[0].append(obj.temperature[0][j])
            temp_k[0].append(obj.k[0][j])
    else:
        for j in range(obj.size[1]):
            x[0].append([obj.boundaries[0], obj.boundaries[0]])
            temp_k[0].append(obj.k[0][j])

    # left boundary
    if obj.boundaries[2] == 0:
        for i in range(obj.size[0]+2):
            x[i].insert(0, x[i][0])
            temp_k[i].insert(0, temp_k[i][0])
    else:
        for i in range(obj.size[0]+2):
            x[i].insert(0, [obj.boundaries[2], obj.boundaries[2]])
            temp_k[i].insert(0, temp_k[i][0])

    # right boundary
    if obj.boundaries[3] == 0:
        for i in range(obj.size[0]+2):
            x[i].append(x[i][-1])
            temp_k[i].append(temp_k[i][-1])
    else:
        for i in range(obj.size[0]+2):
            x[i].append([obj.boundaries[3], obj.boundaries[3]])
            temp_k[i].append(temp_k[i][-1])

    final_x = copy.deepcopy(x)

    # computes
    for i in range(1, obj.size[0]+1):
        for j in range(1, obj.size[1]+1):
            # print(i,j)
            alpha = obj.dt / (2 * obj.rho[i-1][j-1] * obj.Cp[i-1][j-1] *
                              obj.dx * obj.dx)

            gamma = obj.dt / (2 * obj.rho[i-1][j-1] * obj.Cp[i-1][j-1] *
                              obj.dy * obj.dy)

            beta = obj.dt / (obj.rho[i-1][j-1] * obj.Cp[i-1][j-1])

            t_new = ((1 + beta * obj.Q[i-1][j-1]) * x[i][j][0] +
                     alpha * ((temp_k[i+1][j]+temp_k[i][j])*x[i+1][j][0] -
                              (temp_k[i-1][j] + 2 * temp_k[i][j] +
                               temp_k[i+1][j])*x[i][j][0] +
                              (temp_k[i-1][j]+temp_k[i][j])*x[i-1][j][0]) +
                     gamma * ((temp_k[i][j+1]+temp_k[i][j])*x[i][j+1][0] -
                              (temp_k[i][j-1] + 2 * temp_k[i][j] +
                               temp_k[i][j+1])*x[i][j][0] +
                              (temp_k[i][j-1]+temp_k[i][j])*x[i][j-1][0]) +
                     beta * (obj.Q0[i-1][j-1] - obj.Q[i-1][j-1] *
                             obj.amb_temperature))
            final_x[i][j][1] = t_new

    x = final_x
    # updates temperature for next iteration
    for i in range(len(x)):
        for j in range(len(x[1])):
            x[i][j][0] = x[i][j][1]

    # remove boundaries
    removed_temp = x.pop(0)
    removed_temp = x.pop()
    for i in range(len(x)):
        removed_temp = x[i].pop(0)
        removed_temp = x[i].pop()
    del removed_temp

    nx = []
    for i in range(obj.size[0]):
        nx.append([])
        for j in range(obj.size[1]):
            nx[-1].append(x[i][j][0])

    # print(x)
    # print('\n')
    # latent heat
    lheat = copy.copy(obj.lheat)
    for i in range(obj.size[0]):
        for j in range(obj.size[1]):
            k = 0
            for lh in obj.latent_heat[i][j]:
                temper = obj.temperature[i][j][0]
                value = nx[i][j] > lh[0]
                if value and temper <= lh[0] and lheat[i][j][k][1] != lh[1]:
                    en = (obj.Cp[i][j] * obj.rho[i][j] *
                          (nx[i][j] - obj.temperature[i][j][0]))
                    if en + lheat[i][j][k][1] >= lh[1]:
                        lheat[i][j][k][1] = lh[1]
                        energy_temp = lheat[i][j][k][1] + en - lh[1]
                        nx[i][j] = obj.temperature[i][j][0] + \
                            energy_temp / (obj.Cp[i][j] * obj.rho[i][j])
                    else:
                        lheat[i][j][k][1] += en
                        nx[i][j] = obj.temperature[i][j][0]
                value = nx[i][j] < lh[0]
                if value and temper >= lh[0] and lheat[i][j][k][1] != 0:
                    en = (obj.Cp[i][j] * obj.rho[i][j] *
                          (nx[i][j] - obj.temperature[i][j][0]))
                    if en + lheat[i][j][k][1] <= 0.:
                        lheat[i][j][k][1] = 0.
                        energy_temp = (en + lheat[i][j][k][1])
                        nx[i][j] = obj.temperature[i][j][0] + \
                            energy_temp / (obj.Cp[i][j] * obj.rho[i][j])
                    else:
                        lheat[i][j][k][1] += en
                        nx[i][j] = obj.temperature[i][j][0]
                k += 1

    y = copy.deepcopy(obj.temperature)

    # updates the temperature list
    for i in range(obj.size[0]):
        for j in range(obj.size[1]):
            y[i][j][1] = nx[i][j]
            y[i][j][0] = nx[i][j]

    return y, lheat
