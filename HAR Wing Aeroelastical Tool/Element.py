import numpy as np
import math


class Element:
    def __init__(self, node1, node2, node3, rot, length, Ksection, Csection):
        # rot must be in radians!
        self.strainVec = None
        self.strainPVec = None
        self.h0 = None
        self.expG = None
        self.exp2G = None
        self.dexpGdEps = None
        self.dexp2GdEps = None
        self.Jhepp = None

        self.memberJhep = None
        self.memberJpep = None
        self.memberJthetaep = None
        self.memberB = None
        self.memberJhb = None
        self.memberJpb = None
        self.memberJthetab = None
        self.strainm = None
        self.strainpm = None
        self.strainppm = None
        self.lambdm = None

        self.node1 = node1
        M1 = node1.massMatrix

        self.node2 = node2
        M2 = node2.massMatrix

        self.node3 = node3
        M3 = node3.massMatrix

        self.length = length

        self.K = self.length * Ksection  # TODO rigidit
        self.C = self.length * Csection  # TODO damping

        self.Melem = 0.5 * self.length * np.block(
            [[(1 / 4 * M1 + 1 / 12 * M2), (1 / 12 * M1 + 1 / 12 * M2), np.zeros((12, 12))],
             [(1 / 12 * M1 + 1 / 12 * M2), (1 / 12 * M1 + 1 / 2 * M2 + 1 / 12 * M3), (1 / 12 * M2 + 1 / 12 * M3)],
             [np.zeros((12, 12)), (1 / 12 * M2 + 1 / 12 * M3), (1 / 12 * M2 + 1 / 4 * M3)]])  # eq 3.34
        # TODO verificar se precisa add rigid

        # normal N

        N1 = node1.massMatrix[:, 0:3].copy()
        N2 = node2.massMatrix[:, 0:3].copy()
        N3 = node3.massMatrix[:, 0:3].copy()

        self.Nelem = 0.5 * self.length * np.block(
            [[1 / 3 * N1 + 1 / 6 * N2],
             [1 / 6 * N1 + 2 / 3 * N2 + 1 / 6 * N3],
             [1 / 6 * N2 + 1 / 3 * N3]])  # eq 3.49
        # TODO verificar se precisa add rigid

        geomRot = np.array([[1., 0., 0.],
                            [0., math.cos(rot.twist), -math.sin(rot.twist)],
                            [0., math.sin(rot.twist), math.cos(rot.twist)]]) @ \
                  np.array([[math.cos(rot.sweep), -math.sin(rot.sweep), 0.],
                            [math.sin(rot.sweep), math.cos(rot.sweep), 0.],
                            [0., 0., 1.]]) @ \
                  np.array([[math.cos(rot.dihedral), 0., math.sin(rot.dihedral)],
                            [0., 1., 0.],
                            [-math.sin(rot.dihedral), 0., math.cos(rot.dihedral)]])  # eq A.13
        self.elemRot = np.block([[np.eye(3), np.zeros((3, 9))],
                                 [np.zeros((3, 3)), np.eye(3) * geomRot[0, 0], np.eye(3) * geomRot[0, 1],
                                  np.eye(3) * geomRot[0, 2]],
                                 [np.zeros((3, 3)), np.eye(3) * geomRot[1, 0], np.eye(3) * geomRot[1, 1],
                                  np.eye(3) * geomRot[1, 2]],
                                 [np.zeros((3, 3)), np.eye(3) * geomRot[2, 0], np.eye(3) * geomRot[2, 1],
                                  np.eye(3) * geomRot[2, 2]]])  # eq A.12

    # ------------- # ------------- METHODS ------------- # ------------- #

    def setStrain(self, strainVec, strainPVec):
        self.strainVec = strainVec
        self.strainPVec = strainPVec  # TODO descobrir e comentar

    def setH0(self, h0):  # h0 is the vector of initial deformations/deflexions
        self.h0 = self.elemRot @ h0  # eq A.1


