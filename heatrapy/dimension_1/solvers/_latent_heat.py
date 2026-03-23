"""Shared latent heat computation for 1D solvers."""

import copy


def apply_latent_heat(nx, obj):
    """Apply latent heat corrections to computed temperatures.

    Handles phase transitions by absorbing/releasing energy when the
    temperature crosses a transition threshold.

    Parameters
    ----------
    nx : list
        New temperatures for each grid point (modified in place).
    obj : Object
        Thermal object with current state.

    Returns
    -------
    lheat : list
        Updated latent heat accumulation state.

    """
    lheat = copy.copy(obj.lheat)
    for i in range(1, obj.num_points - 1):
        j = 0
        for lh in obj.latent_heat[i]:
            temper = obj.temperature[i][0]
            # heating: crossing transition from below
            if nx[i] > lh[0] and temper <= lh[0] and lheat[i][j][1] != lh[1]:
                en = obj.Cp[i] * obj.rho[i] * (nx[i] - temper)
                if en + lheat[i][j][1] >= lh[1]:
                    lheat[i][j][1] = lh[1]
                    energy_temp = lheat[i][j][1] + en - lh[1]
                    nx[i] = temper + energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    nx[i] = temper
            # cooling: crossing transition from above
            if nx[i] < lh[0] and temper >= lh[0] and lheat[i][j][1] != 0:
                en = obj.Cp[i] * obj.rho[i] * (nx[i] - temper)
                if en + lheat[i][j][1] <= 0.:
                    lheat[i][j][1] = 0.
                    energy_temp = (en + lheat[i][j][1])
                    nx[i] = temper + energy_temp / (obj.Cp[i] * obj.rho[i])
                else:
                    lheat[i][j][1] += en
                    nx[i] = temper
            j += 1

    return lheat
