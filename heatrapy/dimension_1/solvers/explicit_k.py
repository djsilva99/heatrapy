"""Contains the 1D explicit_k solver.

Used to compute unidimensional thermal processes.

"""

import numpy as np

from ._latent_heat import apply_latent_heat


def explicit_k(obj):
    """explicit_k solver.

    Used to compute one time step of 1D systems with x-dependent thermal
    conductivity.

    """
    n = obj.num_points
    s = slice(1, n - 1)

    # extract current temperatures and material properties as arrays
    T = np.array([obj.temperature[i][0] for i in range(n)], dtype=float)
    k = np.array([obj.k[i] if obj.k[i] is not None else 0.0 for i in range(n)])
    rho = np.array([obj.rho[i] if obj.rho[i] is not None else 1.0
                     for i in range(n)])
    Cp = np.array([obj.Cp[i] if obj.Cp[i] is not None else 1.0
                    for i in range(n)])
    Q = np.array([obj.Q[i] if obj.Q[i] is not None else 0.0
                   for i in range(n)])
    Q0 = np.array([obj.Q0[i] if obj.Q0[i] is not None else 0.0
                    for i in range(n)])

    # vectorized FDM stencil for interior points (k varies with x)
    eta = obj.dt / (2.0 * rho[s] * Cp[s] * obj.dx * obj.dx)
    beta = obj.dt / (rho[s] * Cp[s])

    nx = T.copy()
    nx[s] = ((1 + beta * Q[s]) * T[s] +
             eta * ((k[2:n] + k[s]) * T[2:n] -
                    (k[0:n-2] + k[2:n] + 2 * k[s]) * T[s] +
                    (k[0:n-2] + k[s]) * T[0:n-2]) +
             beta * (Q0[s] - Q[s] * obj.amb_temperature))

    # boundaries
    if obj.boundaries[0] == 0:
        nx[0] = T[1]
    else:
        nx[0] = obj.boundaries[0]

    if obj.boundaries[1] == 0:
        nx[n - 1] = T[n - 2]
    else:
        nx[n - 1] = obj.boundaries[1]

    # latent heat (per-element, branching logic)
    nx_list = nx.tolist()
    lheat = apply_latent_heat(nx_list, obj)

    # pack into [current, next] pairs expected by the caller
    y = [[nx_list[i], nx_list[i]] for i in range(n)]

    return y, lheat
