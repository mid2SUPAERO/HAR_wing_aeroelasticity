import numpy as np
import math


def makeDiag(matrix):
    X = np.block([[matrix, np.zeros((4, 4)), np.zeros((4, 4))],
                  [np.zeros((4, 4)), matrix, np.zeros((4, 4))],
                  [np.zeros((4, 4)), np.zeros((4, 4)), matrix]])

    return X


def getKEps(strain):  # schearer pag 28
    ex = strain[0, 0]
    kx = strain[0, 1]
    ky = strain[0, 2]
    kz = strain[0, 3]
    exb = (1 + ex)
    KEps = np.array([[0, exb, 0, 0],
                     [0, 0, kz, -ky],
                     [0, -kz, 0, kx],
                     [0, ky, -kx, 0]])

    return KEps


def getKEps2(strain):
    ex = strain[0, 0]
    kx = strain[0, 1]
    ky = strain[0, 2]
    kz = strain[0, 3]
    exb = (1 + ex)
    lamb2 = (kx ** 2 + ky ** 2 + kz ** 2)
    KEps2 = np.array([[0, 0, exb * kz, -exb * ky],
                      [0, -(lamb2 - kx ** 2), ky * kx, kz * kx],
                      [0, ky * kx, -(lamb2 - ky ** 2), kz * ky],
                      [0, kz * kx, kz * ky, -(lamb2 - kz ** 2)]])

    return KEps2


def getKEps3(strain):
    ex = strain[0, 0]
    kx = strain[0, 1]
    ky = strain[0, 2]
    kz = strain[0, 3]
    exb = (1 + ex)
    lamb2 = (kx ** 2 + ky ** 2 + kz ** 2)
    KEps3 = np.array([[0, -exb * (kz ** 2 + ky ** 2), exb * ky * kx, exb * kz * kx],
                      [0, 0, -kz * lamb2, ky * lamb2],
                      [0, kz * lamb2, 0, -kx * lamb2],
                      [0, -ky * lamb2, kx * lamb2, 0]])

    return KEps3


def getExpG(strain, length):  # eq A8
    KEps = getKEps(strain)
    KEps2 = getKEps2(strain)
    KEps3 = getKEps3(strain)
    ex = strain[0, 0]
    kx = strain[0, 1]
    ky = strain[0, 2]
    kz = strain[0, 3]

    lamb2 = kx ** 2 + ky ** 2 + kz ** 2
    lamb = math.sqrt(lamb2)

    if lamb2 > 1e-10:
        a = (length / lamb2 - math.sin(lamb * length) / lamb ** 3)
        b = (1 - math.cos(lamb * length)) / lamb2
        c = length

    else:
        a = length ** 3 / 6
        b = length ** 2 / 2
        c = length

    expGast = a * KEps3 + b * KEps2 + c * KEps + np.eye(4)
    Th = np.zeros((12, 12))
    Th[0, 0] = 1
    Th[1, 4] = 1
    Th[2, 8] = 1
    Th[3, 1] = 1
    Th[4, 5] = 1
    Th[5, 9] = 1
    Th[6, 2] = 1
    Th[7, 6] = 1
    Th[8, 10] = 1
    Th[9, 3] = 1
    Th[10, 7] = 1
    Th[11, 11] = 1
    expG = Th @ makeDiag(expGast) @ Th.transpose()

    if lamb2 > 1e-10:
        dbdEps = (lamb * length * math.sin(lamb * length) - 2 * (1 - math.cos(lamb * length))) / lamb ** 4 * np.block(
            [[0 * np.eye(4)], [kx * np.eye(4)], [ky * np.eye(4)], [kz * np.eye(4)]])

        dadEps = (3 * math.sin(lamb * length) / lamb ** 5 - (
                2 * length + length * math.cos(lamb * length)) / lamb ** 4) * np.block(
            [[0 * np.eye(4)], [kx * np.eye(4)], [ky * np.eye(4)], [kz * np.eye(4)]])

    else:
        dbdEps = np.zeros((16, 4))
        dadEps = np.zeros((16, 4))

    dKEpsdex = np.zeros((4, 4))
    dKEpsdex[0, 1] = 1
    dKEpsdkx = np.zeros((4, 4))
    dKEpsdkx[2, 3] = 1
    dKEpsdkx[3, 2] = -1
    dKEpsdky = np.zeros((4, 4))
    dKEpsdky[3, 1] = 1
    dKEpsdky[1, 3] = -1
    dKEpsdkz = np.zeros((4, 4))
    dKEpsdkz[1, 2] = 1
    dKEpsdkz[2, 1] = -1  # aparentemente erraram aqui!!! utilizando dKEpsdkz(2,1) = 0, os resultados de frequencia
    # natural ficam exatamente iguais aos das teses

    dKEpsdEps = np.block([[dKEpsdex], [dKEpsdkx], [dKEpsdky], [dKEpsdkz]])
    KEpsident = np.block([[KEps, np.zeros((4, 4)), np.zeros((4, 4)), np.zeros((4, 4))],
                          [np.zeros((4, 4)), KEps, np.zeros((4, 4)), np.zeros((4, 4))],
                          [np.zeros((4, 4)), np.zeros((4, 4)), KEps, np.zeros((4, 4))],
                          [np.zeros((4, 4)), np.zeros((4, 4)), np.zeros((4, 4)), KEps]])

    KEps2ident = np.block([[KEps2, np.zeros((4, 4)), np.zeros((4, 4)), np.zeros((4, 4))],
                           [np.zeros((4, 4)), KEps2, np.zeros((4, 4)), np.zeros((4, 4))],
                           [np.zeros((4, 4)), np.zeros((4, 4)), KEps2, np.zeros((4, 4))],
                           [np.zeros((4, 4)), np.zeros((4, 4)), np.zeros((4, 4)), KEps2]])

    dexpGdEps = dadEps @ KEps3 + a * (
            dKEpsdEps @ KEps2 + KEpsident @ dKEpsdEps @ KEps + KEps2ident @ dKEpsdEps) + dbdEps @ KEps2 + b * (
                        dKEpsdEps @ KEps + KEpsident @ dKEpsdEps) + c * dKEpsdEps

    dexpGdEps = np.block([Th @ makeDiag(dexpGdEps[0:4, 0:4]) @ Th.transpose(),
                          Th @ makeDiag(dexpGdEps[4:8, 0:4]) @ Th.transpose(),
                          Th @ makeDiag(dexpGdEps[8: 12, 0: 4]) @ Th.transpose(),
                          Th @ makeDiag(dexpGdEps[12:16, 0:4]) @ Th.transpose()])

    return expG, dexpGdEps
