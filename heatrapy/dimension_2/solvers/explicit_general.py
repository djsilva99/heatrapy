"""Contains the explicit_general solver.

Used to compute thermal processes

"""

import copy


def explicit_general(obj):
    """explicit_general solver.

    Used to compute one time step of systems with fixed thermal contuctivity.

    """
    x = copy.deepcopy(obj.temperature)

    # bottom boundary
    x.append([])
    if obj.boundaries[3] == 0:
        for j in range(obj.size[1]):
            x[-1].append(obj.temperature[-1][j])
    else:
        # x.append([])
        for j in range(obj.size[1]):
            x[-1].append([obj.boundaries[3], obj.boundaries[3]])

    # top boundary
    x.insert(0, [])
    if obj.boundaries[2] == 0:
        for j in range(obj.size[1]):
            x[0].append(obj.temperature[0][j])
    else:
        for j in range(obj.size[1]):
            x[0].append([obj.boundaries[2], obj.boundaries[2]])

    # left boundary
    if obj.boundaries[0] == 0:
        for i in range(obj.size[0]+2):
            x[i].insert(0, x[i][0])
    else:
        for i in range(obj.size[0]+2):
            x[i].insert(0, [obj.boundaries[0], obj.boundaries[0]])

    # print(len(obj.temperature))
    # right boundary
    if obj.boundaries[1] == 0:
        for i in range(obj.size[0]+2):
            x[i].append(x[i][-1])
    else:
        for i in range(obj.size[0]+2):
            x[i].append([obj.boundaries[1], obj.boundaries[1]])

    # print(obj.temperature)
    # print(x)
    # print(len(x))
    # print(len(obj.temperature))

    final_x = copy.deepcopy(x)
    # print(final_x)
    # print('\n')

    # computes
    for i in range(1, obj.size[0]+1):
        for j in range(1, obj.size[1]+1):
            # print(i, j)
            # print(x[i][j][0])
            # print(x[i][j+1][0])
            # print(x[i][j-1][0])
            # print(x[i-1][j][0])
            # print(x[i+1][j][0])
            # print(obj.Q[i-1][j-1])
            # print(obj.Q)
            # print(obj.k[i-1][j-1])
            # print(obj.rho[i-1][j-1])
            # print(obj.Cp[i-1][j-1])
            # print('\n')
            alpha = obj.dt * \
                obj.k[i-1][j-1] / (obj.rho[i-1][j-1] * obj.Cp[i-1][j-1] *
                               obj.dx * obj.dx)

            gamma = obj.dt * \
                obj.k[i-1][j-1] / (obj.rho[i-1][j-1] * obj.Cp[i-1][j-1] *
                               obj.dy * obj.dy)

            beta = obj.dt / (obj.rho[i-1][j-1] * obj.Cp[i-1][j-1])

            t_new = ((1 + beta * obj.Q[i-1][j-1]) * x[i][j][0] +
                     alpha * (x[i - 1][j][0] - 2 *
                              x[i][j][0] + x[i + 1][j][0]) +
                     gamma * (x[i][j - 1][0] - 2 *
                              x[i][j][0] + x[i][j + 1][0]) +
                     beta * (obj.Q0[i-1][j-1] - obj.Q[i-1][j-1] * obj.amb_temperature))
            final_x[i][j][1] = t_new
            # print(alpha,gamma,beta)
            # print(t_new)
            # print('\n')

    # print(final_x)
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
                if nx[i][j] > lh[0] and temper <= lh[0] and lheat[i][j][k][1] != lh[1]:
                    en = obj.Cp[i][j] * obj.rho[i][j] * (nx[i][j] - obj.temperature[i][j][0])
                    if en + lheat[i][j][k][1] >= lh[1]:
                        lheat[i][j][k][1] = lh[1]
                        energy_temp = lheat[i][j][k][1] + en - lh[1]
                        nx[i][j] = obj.temperature[i][j][0] + \
                            energy_temp / (obj.Cp[i][j] * obj.rho[i][j])
                    else:
                        lheat[i][j][k][1] += en
                        nx[i][j] = obj.temperature[i][j][0]
                if nx[i][j] < lh[0] and temper >= lh[0] and lheat[i][j][k][1] != 0:
                    en = obj.Cp[i][j] * obj.rho[i][j] * (nx[i][j] - obj.temperature[i][j][0])
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
    # print(y)

    # updates the temperature list
    for i in range(obj.size[0]):
        for j in range(obj.size[1]):
            y[i][j][1] = nx[i][j]
            y[i][j][0] = nx[i][j]

    # print(nx)
    return y, lheat










# def explicit_general(obj):
#     """explicit_general solver.

#     Used to compute one time step of systems with fixed thermal contuctivity.

#     """
#     x = copy.deepcopy(obj.temperature)

#     # computes
#     for i in range(1, obj.num_points - 1):

#         alpha = obj.dt * \
#             obj.k[i] / (obj.rho[i] * obj.Cp[i] *
#                         obj.dx * obj.dx)
#         beta = obj.dt / (obj.rho[i] * obj.Cp[i])

#         t_new = ((1 + beta * obj.Q[i]) * obj.temperature[i][0] +
#                  alpha * (obj.temperature[i - 1][0] - 2 *
#                           obj.temperature[i][0] + obj.temperature[i + 1][0]) +
#                  beta * (obj.Q0[i] - obj.Q[i] * obj.amb_temperature))
#         x[i][1] = t_new

#     # left boundary for next time step
#     if obj.boundaries[0] == 0:
#         x[0][1] = obj.temperature[1][1]
#     else:
#         x[0][1] = obj.boundaries[0]

#     # right boundary for next time step
#     if obj.boundaries[1] == 0:
#         x[obj.num_points - 1][1] = obj.temperature[obj.num_points - 2][1]
#     else:
#         x[obj.num_points - 1][1] = obj.boundaries[1]

#     # updates temperature for next iteration
#     for i in range(0, obj.num_points):
#         x[i][0] = x[i][1]

#     nx = []
#     for i in range(len(x)):
#         nx.append(x[i][0])

#     # latent heat
#     lheat = copy.copy(obj.lheat)
#     for i in range(1, obj.num_points - 1):
#         j = 0
#         for lh in obj.latent_heat[i]:
#             temper = obj.temperature[i][0]
#             if nx[i] > lh[0] and temper <= lh[0] and lheat[i][j][1] != lh[1]:
#                 en = obj.Cp[i] * obj.rho[i] * (nx[i] - obj.temperature[i][0])
#                 if en + lheat[i][j][1] >= lh[1]:
#                     lheat[i][j][1] = lh[1]
#                     energy_temp = lheat[i][j][1] + en - lh[1]
#                     nx[i] = obj.temperature[i][0] + \
#                         energy_temp / (obj.Cp[i] * obj.rho[i])
#                 else:
#                     lheat[i][j][1] += en
#                     nx[i] = obj.temperature[i][0]
#             if nx[i] < lh[0] and temper >= lh[0] and lheat[i][j][1] != 0:
#                 en = obj.Cp[i] * obj.rho[i] * (nx[i] - obj.temperature[i][0])
#                 if en + lheat[i][j][1] <= 0.:
#                     lheat[i][j][1] = 0.
#                     energy_temp = (en + lheat[i][j][1])
#                     nx[i] = obj.temperature[i][0] + \
#                         energy_temp / (obj.Cp[i] * obj.rho[i])
#                 else:
#                     lheat[i][j][1] += en
#                     nx[i] = obj.temperature[i][0]
#             j += 1

#     y = copy.deepcopy(obj.temperature)

#     # updates the temperature list
#     for i in range(obj.num_points):
#         y[i][1] = nx[i]
#         y[i][0] = nx[i]

#     return y, lheat
