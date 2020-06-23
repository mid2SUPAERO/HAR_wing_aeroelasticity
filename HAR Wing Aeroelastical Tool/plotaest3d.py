import numpy as np
from math import floor
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt


def plotaest3d(member, translate, ax, option = 'surface'):
    if member.elementsVector[1].node3.geometry is None:
        print('WARNING: Can`t plot structure without geometric data')
    else:
        membersize =len(member.elementsVector)
        j = 0
        deflexao = np.zeros((12, 3 * membersize))
        X = np.zeros((3*membersize, 3))
        Y = np.zeros((3*membersize, 3))
        Z = np.zeros((3*membersize, 3))

        for i in range(0, membersize):
            deflexao[:, j] = member.elementsVector[i].node1.h[:, 0]
            j = j + 1

            deflexao[:, j] = member.elementsVector[i].node2.h[:, 0]
            j = j + 1

            deflexao[:, j] = member.elementsVector[i].node3.h[:, 0]
            j = j + 1

        for i in range(0, deflexao.shape[1]):
            if i % 3 == 0:
                a = member.elementsVector[floor(i / 3)].node1.geometry.a
                b = member.elementsVector[floor(i / 3)].node1.geometry.b
            elif i % 3 == 1:
                a = member.elementsVector[floor(i / 3)].node2.geometry.a
                b = member.elementsVector[floor(i / 3)].node2.geometry.b
            else:
                a = member.elementsVector[floor(i / 3)].node3.geometry.a
                b = member.elementsVector[floor(i / 3)].node3.geometry.b

            X[i, 0] = deflexao[0, i] + translate[0]
            X[i, 1] = deflexao[0, i] + translate[0]
            X[i, 2] = deflexao[0, i] + translate[0]

            Y[i, 0] = deflexao[1, i] + deflexao[7, i] * (b * a + b) + translate[1]
            Y[i, 1] = deflexao[1, i] + translate[1]
            Y[i, 2] = deflexao[1, i] + deflexao[7, i] * (b * a - b) + translate[1]

            Z[i, 0] = deflexao[2, i] + deflexao[8, i] * (b * a + b) + translate[2]
            Z[i, 1] = deflexao[2, i] + translate[2]
            Z[i, 2] = deflexao[2, i] + deflexao[8, i] * (b * a - b) + translate[2]

        # ax.set_title('Aircraft wing deformation in flight')

        # ax.set_xlim3d(-20, 20)
        ax.set_ylim3d(-8*np.max(Y), 8*np.max(Y))
        # ax.set_zlim3d(0, 1.5)
        
        ax.set_zlabel('Displacement (m)')

        if option == 'surface':
            surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis')
        elif option == 'wireframe':
            surf = ax.plot_wireframe(X, Y, Z, color='black')

        return surf
