import numpy as np
from math import floor
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt


def plotaest3d(member, translate):
    if member.elementsVector[1].node3.geometry is None:
        print('WARNING: Can`t plot structure without geometric data')
    else:
        membersize = member.elementsVector.shape[1]
        j = 0
        deflexao = np.zeros((9, 3 * membersize))
        X = np.zeros((membersize, 3))
        Y = np.zeros((membersize, 3))
        Z = np.zeros((membersize, 3))

        for i in range(0, membersize):
            deflexao[:, j] = member.elementsVector[i].node1.h
            j = j + 1

            deflexao[:, j] = member.elementsVector[i].node2.h
            j = j + 1

            deflexao[:, j] = member.elementsVector[i].node3.h
            j = j + 1

        for i in range(0, deflexao.shape[1]):
            if i % 3 == 0:
                a = member[floor(i / 3)].node1.geometry.a
                b = member[floor(i / 3)].node1.geometry.b
            elif i % 3 == 1:
                a = member[floor(i / 3)].node2.geometry.a
                b = member[floor(i / 3)].node2.geometry.b
            else:
                a = member[floor(i / 3)].node3.geometry.a
                b = member[floor(i / 3)].node3.geometry.b

            X[i, 0] = deflexao[0, i] + translate[0]
            X[i, 1] = deflexao[0, i] + translate[0]
            X[i, 2] = deflexao[0, i] + translate[0]

            Y[i, 0] = deflexao[1, i] + deflexao[7, i] * (b * a + b) + translate[1]
            Y[i, 1] = deflexao[1, i] + translate[1]
            Y[i, 2] = deflexao[1, i] + deflexao[7, i] * (b * a - b) + translate[1]

            Z[i, 0] = deflexao[2, i] + deflexao[8, i] * (b * a + b) + translate[2]
            Z[i, 1] = deflexao[2, i] + translate[2]
            Z[i, 2] = deflexao[2, i] + deflexao[8, i] * (b * a - b) + translate[2]

        ax = plt.axes(projection='3d')

        ax.set_title('surface')

        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis')

        return surf
