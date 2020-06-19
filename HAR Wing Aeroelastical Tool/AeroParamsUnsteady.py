import numpy as np
import math


class AeroParamsUnsteady:
    def __init__(self, N, b, a, alpha0, clAlpha, clDelta, cm0, cmDelta, cd0, nDelta):
        self.N = N
        self.b = b
        self.a = a
        self.n = 2  # ?
        self.m = 2  # ?

        self.alpha0 = alpha0
        self.clalpha = clAlpha
        self.cldelta = clDelta
        self.cm0 = cm0
        self.cmdelta = cmDelta
        self.cd0 = cd0
        self.ndelta = nDelta

        self.A = np.zeros((self.N, self.N))
        self.B = None
        self.Bm = None
        self.B1 = None
        self.B2 = None
        self.C = np.zeros((self.m, self.m))
        self.T = np.zeros((self.m, self.n))
        self.K = np.zeros((self.m, self.m))

        AeroParamsUnsteady.calcMatrix(self)

        self.E1 = -np.linalg.inv(self.A)
        self.E2 = np.linalg.inv(self.A) @ self.B1[:, 0]
        self.E3 = np.linalg.inv(self.A) @ self.B1[:, 1]
        self.E4 = np.linalg.inv(self.A) @ self.B2[:, 1]

    # ------------- # ------------- METHODS ------------- # ------------- #
    # TODO criar setmembermatrix em members, esta em Peter.m necessario?
    def calcMatrix(self):  # eq 2.63

        self.T[0, 0] = self.b

        if self.n >= 2:  # TODO verificar se não existe so essa condição
            self.T[0, 1] = -self.b * self.a
            self.T[1, 1] = self.b

        # TODO necessario? if self.n>=3:

        for i1 in range(0, self.m):
            self.K[0, i1] = i1
            self.K[i1, i1] = -i1 / 2

        if self.m >= 1:
            self.C[0, 0] = 1.
        if self.m >= 2:
            self.C[0, 1] = 1.
            self.C[1, 0] = -1. / 2
        # TODO necessario? if m>=3?

        c2 = np.zeros((self.N, 1))

        for i1 in range(0, self.N):
            c2[i1] = 2. / (i1+1)

        s2 = np.zeros((self.m, 1))
        s3 = s2

        for i1 in range(0, self.m):
            if i1 == 0:  # TODO isso faz sentido?
                s2[i1] = 1
            elif i1 == 1:
                s2[i1] = 1. / 2

            s3[i1] = i1

        self.B1 = c2 @ s2.transpose() @ self.T
        self.B2 = c2 @ s3.transpose() @ self.T

        self.A = np.zeros((self.N, self.N))

        D = np.zeros((self.N, self.N))
        b = np.zeros((self.N, 1))
        c = np.zeros((self.N, 1))
        d = np.zeros((self.N, 1))

        for i1 in range(0, self.N):
            for i2 in range(0, self.N):
                if i1 == i2 + 1:
                    D[i1, i2] = 1 / (2 * (i1+1))
                elif i1 == i2 - 1:
                    D[i1, i2] = -1 / (2 * (i1+1))

            if i1 != self.N-1:
                b[i1] = ((-1) ** i1) * (math.factorial(self.N + i1) / math.factorial(self.N - 2 - i1)) * (
                        1 / math.factorial(i1+1) ** 2)

            if i1 == self.N-1:  # TODO else?
                b[i1] = (-1) ** (self.N + 1)

            c[i1] = 2 / (i1 + 1)

            if i1 == 0:  # TODO faz sentido? necessario?
                d[i1] = 1. / 2

        self.B = b
        self.Bm = 1 / 2 * self.C @ np.block([[1], [np.zeros((self.m - 1, 1))]]) @ b.transpose()  # TODO esta certo?
        self.A = D + d @ b.transpose() + c @ d.transpose() + (1 / 2) * c @ b.transpose()

