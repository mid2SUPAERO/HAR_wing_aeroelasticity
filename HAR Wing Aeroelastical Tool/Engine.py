import numpy as np
import math


class Engine:
    def __init__(self, numPI, numMEMB, numELEM, numNODE, Fmax = 0, rho0 = 1, V0 = 1, nrho = 0, nv = 0, alphaf = 0, betaf = 0):
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

# TODO falta corrigir essa função

'''
    def getFPROP(self, ap, manete, rho, U): 
        FPROP = np.zeros((ap.NUMele * 9, 1))



        for i = 1:size(prop, 2)
            Cwb = np.array([[1, 0, 0], [0, math.cos(self.alphaf), -math.sin(self.alphaf)],[ 0, math.sin(self.alphaf), math.cos(self.alphaf)]])@np.array([[math.cos(self.betaf), -math.sin(self.betaf), 0],[ math.sin(self.betaf), math.cos(self.betaf), 0],[ 0, 0, 1]])
            CBb = np.eye(3) @ Cwb

        FPROP[(3 * (prop[i].NODEpos)): (3 + 3 * (prop[i].NODEpos))] = CBb * np.array([[0],[manete(prop(i).numPI) * Fmax * (U / V0) ^ nv * (rho / rho0) ^ nrho],[0]])
'''
