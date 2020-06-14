import numpy as np


class Fuselage:
    def __init__(self, m, pcm, I, ttpcm):
        self.m = m
        self.pcm = pcm
        self.I = I
        self.ttpcm = matrixcross(pcm)  # TODO corrigir
        self.tpcm = self.ttpcm.transpose()
        self.MRB = np.block([[np.eye(3) * m, m * self.ttpcm], [m * self.tpcm, self.I]])  # TODO verificar se não é @
        self.N = self.MRB[:, 0: 3]

    def getCRB(self, omega):
        CRB = np.block([[self.m * omega, self.m * omega * self.ttpcm], [self.m * self.tpcm * omega, self.I * omega]])

        return CRB
