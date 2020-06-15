import numpy as np


def OMEGA(omega, nodesnumber):
    saida = np.zeros((4 * nodesnumber, 4 * nodesnumber))
    for i in range(0, nodesnumber * 4):
        saida[(i * 3): (3 + i * 3), (i * 3): (3 + i * 3)] = omega

    return saida
