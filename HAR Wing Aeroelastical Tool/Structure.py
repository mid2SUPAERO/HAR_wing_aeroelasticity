import math
import numpy as np
from scipy.linalg import expm
from matrixcross import matrixcross

from Node import Node
from Element import Element
from AeroParamsUnsteady import AeroParamsUnsteady
from getExpG import getExpG


def veclin(i):
    vec = np.array([0, 0, 0, 0])
    vec[i] = 1

    return vec


def hdiag(h0):  # criar lambda?
    matrix = np.block([[h0, np.zeros((12, 3))],
                       [np.zeros((12, 1)), h0, np.zeros((12, 2))],
                       [np.zeros((12, 2)), h0, np.zeros((12, 1))],
                       [np.zeros((12, 3)), h0]])

    return matrix


def getJhepep(strainVec, h0, length):
    Jhepep = np.zeros((12, 4, 4))
    delta = 1e-5

    for i in range(0, 4):
        JHepp = getJhep(strainVec + delta * veclin(i), h0, length)
        JHepm = getJhep(strainVec - delta * veclin(i), h0, length)
        Jhepep[:, :, i] = (JHepp - JHepm) / (2 * delta)

    return Jhepep


def getJhep(strainVec, h0, length):  # TODO mudar nome
    expG, dexpGdEps = getExpG(strainVec, length / 2)
    Jhep = dexpGdEps @ np.block([[h0, np.zeros((12, 3))],
                                 [np.zeros((12, 1)), h0, np.zeros((12, 2))],
                                 [np.zeros((12, 2)), h0, np.zeros((12, 1))],
                                 [np.zeros((12, 3)), h0]])

    return Jhep


def getnumJhep(strainVec, h0, length):  # TODO mudar nome
    delta = 1e-5
    Jhep = np.zeros((1000, 4))  # TODO rodar testes e corrigir o valor de linhas IMPORTANTE!!!

    for i in range(0, 4):
        expG, dexpGdeps = getnumexpG(strainVec + delta * veclin(i), length / 2)
        expGm, dexpGdepsm = getnumexpG(strainVec - delta * veclin(i), length / 2)
        exp2G, dexp2Gdeps = getnumexpG(strainVec + delta * veclin(i), length)
        exp2Gm, dexp2Gdepsm = getnumexpG(strainVec - delta * veclin(i), length)
        h0p = h0
        h0m = h0
        h1p = expG @ h0
        h1m = expGm @ h0
        h2p = exp2G @ h0
        h2m = exp2Gm @ h0
        Jhep[:, i] = (np.block([[h0p], [h1p], [h2p]]) - np.block([[h0m], [h1m], [h2m]])) / (2 * delta)

    return Jhep


def getnumexpG(strainVec, length):  # TODO mudar nome
    ex = strainVec[0]
    kx = strainVec[1] * np.eye(3)
    ky = strainVec[2] * np.eye(3)
    kz = strainVec[3] * np.eye(3)
    exb = (1 + ex) * np.eye(3)
    KEps = np.block([[np.zeros((3, 3)), exb, np.zeros((3, 3)), np.zeros((3, 3))],
                     [np.zeros((3, 3)), np.zeros((3, 3)), kz, -ky],
                     [np.zeros((3, 3)), -kz, np.zeros((3, 3)), kx],
                     [np.zeros((3, 3)), ky, -kx, np.zeros((3, 3))]])
    expG = expm(KEps * length)
    dexpGdeps = []

    return expG, dexpGdeps


def product3dXvec(A, v):  # TODO checar IMPORTANTE!!!
    matrix = np.zeros((A.shape[1], A.shape[2]))

    for i in range(0, A.shape[0]):
        matrix[:, i] = np.reshape(A[:, i, :], (36, 4), order="F") @ v.tranpose()

    return matrix


class Structure:
    def __init__(self, numOfElems):
        self.numOfElems = numOfElems
        self.elemPos = 0
        self.elementsVector = [0] * self.numOfElems  # np.zeros((self.numOfElems, 1))
        self.Jhep = None
        self.Jpep = None
        self.Jthetaep = None
        self.B = None
        self.Jhb = None
        self.Jpb = None
        self.Jthetab = None

    # ------------- # ------------- METHODS ------------- # ------------- #

    def addElem(self, elem):
        self.elementsVector[self.elemPos] = elem
        self.elementsVector[self.elemPos].setStrain(np.zeros((1, 4)), np.zeros((1, 4)))
        self.elemPos += 1  # TODO adicionar erro se elemenPos>numOfElems

    def update(self):
        Structure.setExponentials(self)
        Structure.setHNodes(self)

    def setExponentials(self):
        for ii in range(0, self.numOfElems):
            self.elementsVector[ii].expG, self.elementsVector[ii].dexpGdEps = getExpG(self.elementsVector[ii].strainVec,
                                                                                      self.elementsVector[
                                                                                          ii].length / 2)
            self.elementsVector[ii].exp2G, self.elementsVector[ii].dexp2GdEps = getExpG(
                self.elementsVector[ii].strainVec,
                self.elementsVector[ii].length)

    def setHNodes(self):
        for ii in range(0, self.numOfElems):
            if ii != 0:
                self.elementsVector[ii].h0 = self.elementsVector[ii].elemRot @ self.elementsVector[ii-1].node3.h

            self.elementsVector[ii].node1.h = self.elementsVector[ii].h0
            self.elementsVector[ii].node2.h = self.elementsVector[ii].expG @ self.elementsVector[ii].h0
            self.elementsVector[ii].node3.h = self.elementsVector[ii].exp2G @ self.elementsVector[ii].h0

    def setJhep2(self):  # TODO mudar nome
        for ii in range(0, self.numOfElems):
            h0 = hdiag(self.elementsVector[ii].h0)
            self.elementsVector[ii].node1.Jhep = np.zeros((12, 4))
            self.elementsVector[ii].node2.Jhep = self.elementsVector[ii].dexpGdEps @ h0
            self.elementsVector[ii].node3.Jhep = self.elementsVector[ii].dexp2GdEps @ h0
            self.elementsVector[ii].Jhep = np.block([[self.elementsVector[ii].node1.Jhep],
                                                     [self.elementsVector[ii].node2.Jhep],
                                                     [self.elementsVector[ii].node3.Jhep]])

    def setJhepp(self):
        for ii in range(0, self.numOfElems):
            Jhepep1 = getJhepep(self.elementsVector[ii].h0, 0)
            Jhepep2 = getJhepep(self.elementsVector[ii].h0,
                                self.elementsVector[ii].length)
            Jhepep3 = getJhepep(self.elementsVector[ii].h0,
                                2 * self.elementsVector[ii].length)
            self.elementsVector[ii].Jhepp = product3dXvec(self.elementsVector[ii].strainPVec)

    def getmemberJhep(self):
        Jhep = np.zeros((self.numOfElems * 12 * 3, self.numOfElems * 4))
        for ii in range(0, self.numOfElems):
            for jj in range(0, self.numOfElems):
                if ii == jj:
                    h0 = hdiag(self.elementsVector[ii].h0)
                    Jhep1 = np.zeros((12, 4))
                    Jhep2 = self.elementsVector[ii].dexpGdEps @ h0
                    Jhep3 = self.elementsVector[ii].dexp2GdEps @ h0
                    Jhep[(ii * 36): (ii * 36 + 36), (jj * 4): (jj * 4 + 4)] = np.block([[Jhep1], [Jhep2], [Jhep3]])

                elif ii > jj:
                    Jhep[(ii * 36): (ii * 36 + 12), (jj * 4): (jj * 4 + 4)] = self.elementsVector[ii].elemRot @ \
                                                                              Jhep[((ii-1) * 36 + 24): ((ii-1) * 36 + 36),
                                                                              (jj * 4): (jj * 4 + 4)]
                    Jhep[(ii * 36 + 12): (ii * 36 + 24), (jj * 4): (jj * 4 + 4)] = self.elementsVector[ii].expG @ \
                                                                                   Jhep[(ii * 36): (ii * 36 + 12),
                                                                                   (jj * 4): (jj * 4 + 4)]
                    Jhep[(ii * 36 + 24): (ii * 36 + 36), (jj * 4): (jj * 4 + 4)] = self.elementsVector[ii].expG @ \
                                                                                   Jhep[(ii * 36 + 12): (ii * 36 + 24),
                                                                                   (jj * 4): (jj * 4 + 4)]

        return Jhep

    def getJhb(self):

        Jhb = np.zeros((self.numOfElems * 12 * 3, 6))
        for ii in range(0, self.numOfElems):
            # node 1
            h = self.elementsVector[ii].node1.h
            p = h[0:3]
            wx = h[3:6]
            wy = h[6:9]
            wz = h[9:12]

            Jhb[(3 * 12 * ii): (3 + 3 * 12 * ii), 0: 3] = np.eye(3)
            Jhb[(3 * 12 * ii): (12 + 3 * 12 * ii), 3: 6] = np.block([[matrixcross(p)],
                                                                     [matrixcross(wx)],
                                                                     [matrixcross(wy)],
                                                                     [matrixcross(wz)]])

            # node 2
            h = self.elementsVector[ii].node2.h
            p = h[0:3]
            wx = h[3:6]
            wy = h[6:9]
            wz = h[9:12]

            Jhb[(12 + 3 * 12 * ii): (15 + 3 * 12 * ii), 0: 3] = np.eye(3)
            Jhb[(12 + 3 * 12 * ii): (24 + 3 * 12 * ii), 3: 6] = np.block([[matrixcross(p)],
                                                                          [matrixcross(wx)],
                                                                          [matrixcross(wy)],
                                                                          [matrixcross(wz)]])

            # node 3
            h = self.elementsVector[ii].node3.h
            p = h[0:3]
            wx = h[3:6]
            wy = h[6:9]
            wz = h[9:12]

            Jhb[(24 + 3 * 12 * ii): (27 + 3 * 12 * ii), 0: 3] = np.eye(3)
            Jhb[(24 + 3 * 12 * ii): (36 + 3 * 12 * ii), 3: 6] = np.block([[matrixcross(p)],
                                                                          [matrixcross(wx)],
                                                                          [matrixcross(wy)],
                                                                          [matrixcross(wz)]])

        return Jhb

    def getJpb(self, Jhb):
        Jpbeta = np.zeros((self.numOfElems * 9, 6))
        for i in range(0, self.numOfElems * 3):
            Jpbeta[(3 * i): (3 + 3 * i), :] = Jhb[(i * 12): (3 + i * 12), :]

        return Jpbeta

    def getJthetab(self):
        Jthetabeta = np.zeros((self.numOfElems * 9, 6))

        for i in range(0, self.numOfElems):
            Jthetabeta[(9 * i): (9 + 9 * i), 3: 6] = np.block([[np.eye(3)], [np.eye(3)], [np.eye(3)]])

        return Jthetabeta

    def getJpep(self, Jhep):
        size = int(Jhep.shape[1] / 4)
        Jpep = np.zeros((size * 3*3, size * 4))

        for i in range(0, size * 3):
            Jpep[(i * 3): (3 + i * 3), :] = Jhep[(i * 12): (3 + i * 12), :]

        return Jpep

    def getMe(self):
        Me = np.zeros((self.numOfElems * 36, self.numOfElems * 36))

        for i in range(0, self.numOfElems):
            Me[(i * 36): (i * 36 + 36), (i * 36): (i * 36 + 36)] = self.elementsVector[i].Melem

        return Me

    def getN(self):
        N = np.zeros((self.numOfElems * 36, 3))

        for i in range(0, self.numOfElems):
            N[(i * 36): (i * 36 + 36), 0: 3] = self.elementsVector[i].Nelem

        return N

    def getK(self):
        K = np.zeros((self.numOfElems * 4, self.numOfElems * 4))

        for i in range(0, self.numOfElems):
            K[(i * 4): (i * 4 + 4), (i * 4): (i * 4 + 4)] = self.elementsVector[i].K

        return K

    def getC(self):
        C = np.zeros((self.numOfElems * 4, self.numOfElems * 4))

        for i in range(0, self.numOfElems):
            C[(i * 4): (i * 4 + 4), (i * 4): (i * 4 + 4)] = self.elementsVector[i].C
        return C

    def getB(self):
        Bfe = np.block([[np.eye(3) * 1 / 3, np.eye(3) * 1 / 6, np.zeros((3, 3))],
                        [np.eye(3) * 1 / 6, np.eye(3) * 2 / 3, np.eye(3) * 1 / 6],
                        [np.zeros((3, 3)), np.eye(3) * 1 / 6, np.eye(3) * 1 / 3]])
        Bf = np.zeros((self.numOfElems * 9, self.numOfElems * 9))

        for i in range(0, self.numOfElems):
            Bf[(i * 9): (i * 9 + 9), (i * 9): (i * 9 + 9)] = 0.5 * self.elementsVector[i].length * Bfe

        return Bf

    def getJthetaep(self, Jhep):
        Jthetaep = np.zeros((self.numOfElems * 9, self.numOfElems * 4))

        for i in range(0, self.numOfElems):
            for ii in range(0, self.numOfElems):
                for j in range(0, 4):
                    if i >= ii:
                        # node 1
                        dtdepz = Jhep[(3 + i * 36):(6 + i * 36), j + ii * 4].transpose() @ self.elementsVector[
                                                                                               i].node1.h[6:9]
                        dtdepx = Jhep[(6 + i * 36):(9 + i * 36), j + ii * 4].transpose() @ self.elementsVector[
                                                                                               i].node1.h[9:12]
                        dtdepy = Jhep[(9 + i * 36):(12 + i * 36), j + ii * 4].transpose() @ self.elementsVector[
                                                                                                i].node1.h[3:6]
                        dtdex1 = np.block(
                            [self.elementsVector[i].node1.h[3:6], self.elementsVector[i].node1.h[6: 9],
                             self.elementsVector[i].node1.h[9: 12]]) @ np.block(
                            [[dtdepx], [dtdepy], [dtdepz]])

                        # node 2
                        dtdepz = Jhep[(3 + 12 + i * 36):(6 + 12 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node2.h[
                                 6:9]
                        dtdepx = Jhep[(6 + 12 + i * 36):(9 + 12 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node2.h[
                                 9:12]
                        dtdepy = Jhep[(9 + 12 + i * 36):(12 + 12 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node2.h[
                                 3:6]
                        dtdex2 = np.block(
                            [self.elementsVector[i].node2.h[3:6], self.elementsVector[i].node2.h[6: 9],
                             self.elementsVector[i].node2.h[9: 12]]) @ np.block(
                            [[dtdepx], [dtdepy], [dtdepz]])

                        # node 3
                        dtdepz = Jhep[(3 + 24 + i * 36):(6 + 24 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node3.h[6:9]
                        dtdepx = Jhep[(6 + 24 + i * 36):(9 + 24 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node3.h[9:12]
                        dtdepy = Jhep[(9 + 24 + i * 36):(12 + 24 + i * 36), j + ii * 4].transpose() @ \
                                 self.elementsVector[i].node3.h[3:6]
                        dtdex3 = np.block([self.elementsVector[i].node3.h[3:6], self.elementsVector[i].node3.h[6: 9],
                                           self.elementsVector[i].node3.h[9: 12]]) @ np.block(
                            [[dtdepx], [dtdepy], [dtdepz]])

                        Jthetaep[(i * 9): (i * 9 + 9), (ii * 4 + j)] = np.block([dtdex1.T, dtdex2.T, dtdex3.T])

        return Jthetaep