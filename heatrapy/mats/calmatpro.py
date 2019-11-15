"""Contains the class calmatpro.

Used to access the properties of materials from the database

"""

import numpy as np


class CalMatPro:
    """CalMatPro class.

    Loads and gives physical properties of materials from the database

    """

    def __init__(self, name_tadi, name_tadd, name_cpa, name_cp0, name_k0,
                 name_ka, name_rho0, name_rhoa, name_lheat0, name_lheata):
        """Object initialization.

        loads the physical properties of a materials
        nameTadi is the adiabatic temperature increase file
        nameTadd is the adiabatic temperature decrease file
        nameCp0 is the specific heat file
        nameCpA is the active specific heat file
        namek0 is the thermal conductivity file
        namekA is the active thermal conductivity file
        namerho0 is the density file
        namerhoA is the active density file

        """
        # adiabatic temperature increase
        input = open(name_tadi, 'r')
        s = input.readlines()
        self.xTadi = []
        self.yTadi = []
        for line in s:
            pair_tad = line.split()
            self.xTadi.append(float(pair_tad[0]))
            self.yTadi.append(float(pair_tad[1]))
        input.close()

        # adiabatic temperature decrease
        input = open(name_tadd, 'r')
        s = input.readlines()
        self.xTadd = []
        self.yTadd = []
        for line in s:
            pair_tad = line.split()
            self.xTadd.append(float(pair_tad[0]))
            self.yTadd.append(float(pair_tad[1]))
        input.close()

        # specific heat
        input = open(name_cpa, 'r')
        s = input.readlines()
        self.xCpA = []
        self.yCpA = []
        for line in s:
            pair_tad = line.split()
            self.xCpA.append(float(pair_tad[0]))
            self.yCpA.append(float(pair_tad[1]))
        input.close()

        # active specific heat
        input = open(name_cp0, 'r')
        s = input.readlines()
        self.xCp0 = []
        self.yCp0 = []
        for line in s:
            pair_tad = line.split()
            self.xCp0.append(float(pair_tad[0]))
            self.yCp0.append(float(pair_tad[1]))
        input.close()

        # thermal conductivity
        input = open(name_k0, 'r')
        s = input.readlines()
        self.xk0 = []
        self.yk0 = []
        for line in s:
            pair_tad = line.split()
            self.xk0.append(float(pair_tad[0]))
            self.yk0.append(float(pair_tad[1]))
        input.close()

        # active thermal conductivity
        input = open(name_ka, 'r')
        s = input.readlines()
        self.xkA = []
        self.ykA = []
        for line in s:
            pair_tad = line.split()
            self.xkA.append(float(pair_tad[0]))
            self.ykA.append(float(pair_tad[1]))
        input.close()

        # density
        input = open(name_rho0, 'r')
        s = input.readlines()
        self.xrho0 = []
        self.yrho0 = []
        for line in s:
            pair_tad = line.split()
            self.xrho0.append(float(pair_tad[0]))
            self.yrho0.append(float(pair_tad[1]))
        input.close()

        # active density
        input = open(name_rhoa, 'r')
        s = input.readlines()
        self.xrhoA = []
        self.yrhoA = []
        for line in s:
            pair_tad = line.split()
            self.xrhoA.append(float(pair_tad[0]))
            self.yrhoA.append(float(pair_tad[1]))
        input.close()

        # latent heat
        input = open(name_lheat0, 'r')
        s = input.readlines()
        self.latent_heat0 = []
        for line in s:
            latent_values = line.split()
            self.latent_heat0.append(
                (float(latent_values[0]), float(latent_values[1]))
            )
        input.close()

        # active latent heat
        input = open(name_lheata, 'r')
        s = input.readlines()
        self.latent_heatA = []
        for line in s:
            latent_values = line.split()
            self.latent_heatA.append(
                (float(latent_values[0]),float(latent_values[1]))
            )
        input.close()

    def tadi(self, temperature):
        """Adiabatic temperature increase.

        Gives the adiabatic temperature increase for a given temperature

        """
        return np.interp(temperature, self.xTadi, self.yTadi)

    def tadd(self, temperature):
        """Adiabatic temperature decrease.

        Gives the adiabatic temperature decrease for a given temperature

        """
        return np.interp(temperature, self.xTadd, self.yTadd)

    def cpa(self, temperature):
        """Active specific heat.

        Gives the active specific heat for a given temperature

        """
        return np.interp(temperature, self.xCpA, self.yCpA)

    def cp0(self, temperature):
        """Inactive specific heat.

        Gives the specific heat for a given temperature

        """
        return np.interp(temperature, self.xCp0, self.yCp0)

    def k0(self, temperature):
        """Inactive thermal conductivity.

        Gives the thermal conductivity for a given temperature

        """
        return np.interp(temperature, self.xk0, self.yk0)

    def ka(self, temperature):
        """Active thermal conductivity.

        Gives the active thermal conductivity for a given temperature

        """
        return np.interp(temperature, self.xkA, self.ykA)

    def rho0(self, temperature):
        """Inactive density.

        Gives the density for a given temperature

        """
        return np.interp(temperature, self.xrho0, self.yrho0)

    def rhoa(self, temperature):
        """Active density.

        Gives the active density for a given temperature

        """
        return np.interp(temperature, self.xrhoA, self.yrhoA)

    def lheat0(self):
        """Inactive latent heat.

        Gives the latent heat list

        """
        return self.latent_heat0

    def lheata(self):
        """Active latent heat.

        Gives the active latent heat list

        """
        return self.latent_heatA
