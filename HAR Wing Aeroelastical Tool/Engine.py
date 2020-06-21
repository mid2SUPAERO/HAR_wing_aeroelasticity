import numpy as np
from math import sin, cos


class Engine:
    def __init__(self, numPI, numMEMB, numELEM, numNODE, Fmax=0, rho0=1, V0=1, nrho=0, nv=0, alphaf=0, betaf=0):
        self.numPI = numPI
        self.numMEMB = numMEMB
        self.numELEM = numELEM
        self.numNODE = numNODE
        self.Fmax = Fmax
        self.rho0 = rho0
        self.V0 = V0
        self.nrho = nrho
        self.nv = nv
        self.alphaf = alphaf
        self.betaf = betaf
        self.NODEpos = None

    def getFPROP(self, ap, manete, rho, U):
        FPROP = np.zeros((ap.NUMele * 9, 1))

        for i in range(0, ap.prop.shape[0]):  # TODO checar esse loop
            Cwb = np.array([[1, 0, 0], [0, cos(self.alphaf), -sin(self.alphaf)],
                            [0, sin(self.alphaf), cos(self.alphaf)]]) @ np.array(
                [[cos(self.betaf), -sin(self.betaf), 0], [sin(self.betaf), cos(self.betaf), 0], [0, 0, 1]])
            CBb = np.eye(3) @ Cwb

            FPROP[int(3 * (ap.prop[i].NODEpos-1)): int(3 + 3 * (ap.prop[i].NODEpos-1))] = CBb @ np.array([[0], [
                manete * self.Fmax * (U[0] / self.V0) ** self.nv * (rho / self.rho0) ** self.nrho],
                                                                                            [0]])

        # originally: manete[ap.prop[i].numPI - 1]

        return FPROP
