import numpy as np
from matrixcross import matrixcross


class Fuselage:
    def __init__(self, m, pcm, I):
        self.m = m
        self.pcm = pcm
        self.I = I
        self.ttpcm = matrixcross(pcm)
        self.tpcm = self.ttpcm.transpose()
        self.MRB = np.block([[np.eye(3) * m, m * self.ttpcm], [m * self.tpcm, self.I]])
        self.N = self.MRB[:, 0: 3]

    def getCRB(self, omega):
        CRB = np.block([[self.m * omega, self.m * omega * self.ttpcm], [self.m * self.tpcm * omega, self.I * omega]])

        return CRB
