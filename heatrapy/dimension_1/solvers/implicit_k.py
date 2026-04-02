"""Contains the 1D implicit_k solver.

Used to compute unidimensional thermal processes

"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ._latent_heat import apply_latent_heat

if TYPE_CHECKING:
    from ..objects.object import Object


def implicit_k(
    obj: Object,
) -> tuple[list[list[float]], list[list[list[float]]]]:
    """implicit_k solver.

    Used to compute one time step of 1D systems with k-dependent
    thermal conductivities.

    Parameters
    ----------
    obj : Object
        Thermal object with current state.

    Returns
    -------
    y : list[list[float]]
        Updated temperatures as ``[[T, T], ...]`` pairs.
    lheat : list[list[list[float]]]
        Updated latent heat state.

    Raises
    ------
    TypeError
        If ``obj`` is not a thermal Object.

    """
    if not hasattr(obj, 'num_points'):
        raise TypeError(
            f"obj must be a thermal Object, "
            f"got {type(obj).__name__}"
        )

    n: int = obj.num_points

    # initializes the matrices for the equation system
    a: np.ndarray = np.zeros((n, n))
    b: np.ndarray = np.zeros(n)

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

    # build tridiagonal system and solve
    for i in range(1, n - 1):
        gamma: float = (
            4. * obj.rho[i] * obj.Cp[i]
            * obj.dx * obj.dx / obj.dt
        )

        a[i][i - 1] = obj.k[i - 1] + obj.k[i]
        a[i][i] = -(
            gamma + obj.k[i + 1] + obj.k[i - 1]
            + 2. * obj.k[i]
            - 2 * obj.dt * obj.dt * obj.Q[i]
        )
        a[i][i + 1] = obj.k[i + 1] + obj.k[i]
        b[i] = (
            -(obj.k[i + 1] + obj.k[i])
            * obj.temperature[i + 1][0]
            + (-gamma + obj.k[i + 1] + obj.k[i - 1]
               + 2. * obj.k[i]
               - 2 * obj.dt * obj.dt * obj.Q[i])
            * obj.temperature[i][0]
            - (obj.k[i - 1] + obj.k[i])
            * obj.temperature[i - 1][0]
            - 4. * obj.dx * obj.dx
            * (obj.Q0[i] - obj.Q[i] * obj.amb_temperature)
        )

    x: np.ndarray = np.linalg.solve(a, b)

    # latent heat
    nx_list: list[float] = x.tolist()
    lheat: list[list[list[float]]] = apply_latent_heat(
        nx_list, obj
    )

    # pack into [current, next] pairs expected by the caller
    y: list[list[float]] = [
        [nx_list[i], nx_list[i]] for i in range(n)
    ]

    return y, lheat
