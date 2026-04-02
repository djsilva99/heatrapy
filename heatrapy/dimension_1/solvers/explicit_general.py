"""Contains the 1D explicit_general solver.

Used to compute unidimensional thermal processes

"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ._latent_heat import apply_latent_heat

if TYPE_CHECKING:
    from ..objects.object import Object


def explicit_general(
    obj: Object,
) -> tuple[list[list[float]], list[list[list[float]]]]:
    """explicit_general solver.

    Used to compute one time step of 1D systems with fixed thermal
    conductivity.

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
    s: slice = slice(1, n - 1)

    # extract current temperatures and material properties
    T: np.ndarray = np.array(
        [obj.temperature[i][0] for i in range(n)], dtype=float
    )
    k: np.ndarray = np.array(
        [obj.k[i] if obj.k[i] is not None else 0.0
         for i in range(n)]
    )
    rho: np.ndarray = np.array(
        [obj.rho[i] if obj.rho[i] is not None else 1.0
         for i in range(n)]
    )
    Cp: np.ndarray = np.array(
        [obj.Cp[i] if obj.Cp[i] is not None else 1.0
         for i in range(n)]
    )
    Q: np.ndarray = np.array(
        [obj.Q[i] if obj.Q[i] is not None else 0.0
         for i in range(n)]
    )
    Q0: np.ndarray = np.array(
        [obj.Q0[i] if obj.Q0[i] is not None else 0.0
         for i in range(n)]
    )

    # vectorized FDM stencil for interior points
    alpha: np.ndarray = (
        obj.dt * k[s] / (rho[s] * Cp[s] * obj.dx * obj.dx)
    )
    beta: np.ndarray = obj.dt / (rho[s] * Cp[s])

    nx: np.ndarray = T.copy()
    nx[s] = (
        (1 + beta * Q[s]) * T[s]
        + alpha * (T[0:n-2] - 2 * T[s] + T[2:n])
        + beta * (Q0[s] - Q[s] * obj.amb_temperature)
    )

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
    nx_list: list[float] = nx.tolist()
    lheat: list[list[list[float]]] = apply_latent_heat(
        nx_list, obj
    )

    # pack into [current, next] pairs expected by the caller
    y: list[list[float]] = [
        [nx_list[i], nx_list[i]] for i in range(n)
    ]

    return y, lheat
