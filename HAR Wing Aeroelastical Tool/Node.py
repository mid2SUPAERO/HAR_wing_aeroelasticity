import numpy as np


# TODO add comments
class Node:
    def __init__(self, linDens, coordCG, inertiaMatrix, length, aeroParams, geometry=None):
        self.aeroParams = None
        self.massMatrix = None
        self.h = None  # TODO renomear
        self.Jhep = None

        self.linDens = linDens
        self.coordCG = coordCG
        self.inertiaMatrix = inertiaMatrix
        self.length = length
        self.geometry = geometry  # TODO verificar depois

        self.aeroParams = aeroParams

        Node.calcMassMatrix(self)

    # ------------- # ------------- METHODS ------------- # ------------- #

    def calcMassMatrix(self):
        linDensMatrix = self.linDens * np.eye(3)

        rx = self.coordCG[0] * np.eye(3)
        ry = self.coordCG[1] * np.eye(3)
        rz = self.coordCG[2] * np.eye(3)

        Ixx = self.inertiaMatrix[0, 0] * np.eye(3)
        Ixy = self.inertiaMatrix[0, 1] * np.eye(3)
        Ixz = self.inertiaMatrix[0, 2] * np.eye(3)
        Iyx = self.inertiaMatrix[1, 0] * np.eye(3)
        Iyy = self.inertiaMatrix[1, 1] * np.eye(3)
        Iyz = self.inertiaMatrix[1, 2] * np.eye(3)
        Izx = self.inertiaMatrix[2, 0] * np.eye(3)
        Izy = self.inertiaMatrix[2, 1] * np.eye(3)
        Izz = self.inertiaMatrix[2, 2] * np.eye(3)

        self.massMatrix = np.block([[linDensMatrix, linDensMatrix @ rx, linDensMatrix @ ry, linDensMatrix @ rz],
                                    [linDensMatrix @ rx, 0.5 * (Iyy + Izz - Ixx), Ixy, Ixz],
                                    [linDensMatrix @ ry, Iyx, 0.5 * (Ixx + Izz - Iyy), Iyz],
                                    [linDensMatrix @ rz, Izx, Izy, 0.5 * (Ixx + Iyy - Izz)]])
