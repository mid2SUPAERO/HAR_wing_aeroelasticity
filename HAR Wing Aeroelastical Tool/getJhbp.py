import numpy as np
from matrixcross import matrixcross


def getJhpb(hp, nodesnumber):
    saida = np.zeros((12 * nodesnumber, 6))

    for i in range(0, nodesnumber):
        saida[(i * 12): (3 + i * 12), 3: 6] = matrixcross(hp[(i * 12): (3 + i * 12)])
        saida[(3 + i * 12): (6 + i * 12), 3: 6] = matrixcross(hp[(3 + i * 12): (6 + i * 12)])
        saida[(6 + i * 12): (9 + i * 12), 3: 6] = matrixcross(hp[(6 + i * 12): (9 + i * 12)])
        saida[(9 + i * 12): (12 + i * 12), 3: 6] = matrixcross(hp[(9 + i * 12): (12 + i * 12)])

    return saida
