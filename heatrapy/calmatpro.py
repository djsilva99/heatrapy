# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the class calmatpro.

Used to access the properties of materials from the database

"""

import numpy as np


class calmatpro:

    """Calmatpro class

    Loads and gives physical properties of materials from the database

    """

    def __init__(self, nameTadi, nameTadd, nameCpA, nameCp0, namek0, namekA,
                 namerho0, namerhoA):
        """Initializes the object.

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
        input = open(nameTadi, 'r')
        s = input.readlines()
        self.xTadi = []
        self.yTadi = []
        for line in s:
            pairTad = line.split()
            self.xTadi.append(float(pairTad[0]))
            self.yTadi.append(float(pairTad[1]))
        input.close()

        # adiabatic temperature decrease
        input = open(nameTadd, 'r')
        s = input.readlines()
        self.xTadd = []
        self.yTadd = []
        for line in s:
            pairTad = line.split()
            self.xTadd.append(float(pairTad[0]))
            self.yTadd.append(float(pairTad[1]))
        input.close()

        # specific heat
        input = open(nameCpA, 'r')
        s = input.readlines()
        self.xCpA = []
        self.yCpA = []
        for line in s:
            pairTad = line.split()
            self.xCpA.append(float(pairTad[0]))
            self.yCpA.append(float(pairTad[1]))
        input.close()

        # active specific heat
        input = open(nameCp0, 'r')
        s = input.readlines()
        self.xCp0 = []
        self.yCp0 = []
        for line in s:
            pairTad = line.split()
            self.xCp0.append(float(pairTad[0]))
            self.yCp0.append(float(pairTad[1]))
        input.close()

        # thermal conductivity
        input = open(namek0, 'r')
        s = input.readlines()
        self.xk0 = []
        self.yk0 = []
        for line in s:
            pairTad = line.split()
            self.xk0.append(float(pairTad[0]))
            self.yk0.append(float(pairTad[1]))
        input.close()

        # active thermal conductivity
        input = open(namekA, 'r')
        s = input.readlines()
        self.xkA = []
        self.ykA = []
        for line in s:
            pairTad = line.split()
            self.xkA.append(float(pairTad[0]))
            self.ykA.append(float(pairTad[1]))
        input.close()

        # density
        input = open(namerho0, 'r')
        s = input.readlines()
        self.xrho0 = []
        self.yrho0 = []
        for line in s:
            pairTad = line.split()
            self.xrho0.append(float(pairTad[0]))
            self.yrho0.append(float(pairTad[1]))
        input.close()

        # active density
        input = open(namerhoA, 'r')
        s = input.readlines()
        self.xrhoA = []
        self.yrhoA = []
        for line in s:
            pairTad = line.split()
            self.xrhoA.append(float(pairTad[0]))
            self.yrhoA.append(float(pairTad[1]))
        input.close()

    def tadi(self, Temperature):
        """gives the adiabatic temperature increase for a given temperature"""

        return np.interp(Temperature, self.xTadi, self.yTadi)

    def tadd(self, Temperature):
        """gives the adiabatic temperature decrease for a given temperature"""

        return np.interp(Temperature, self.xTadd, self.yTadd)

    def cpa(self, Temperature):
        """gives the active specific heat for a given temperature"""

        return np.interp(Temperature, self.xCpA, self.yCpA)

    def cp0(self, Temperature):
        """gives the specific heat for a given temperature"""

        return np.interp(Temperature, self.xCp0, self.yCp0)

    def k0(self, Temperature):
        """gives the thermal conductivity for a given temperature"""

        return np.interp(Temperature, self.xk0, self.yk0)

    def ka(self, Temperature):
        """gives the active thermal conductivity for a given temperature"""

        return np.interp(Temperature, self.xkA, self.ykA)

    def rho0(self, Temperature):
        """method gives the density for a given temperature"""

        return np.interp(Temperature, self.xrho0, self.yrho0)

    def rhoa(self, Temperature):
        """method gives the active density for a given temperature"""

        return np.interp(Temperature, self.xrhoA, self.yrhoA)
