import numpy as np

def matrixcross(p):
    matrix = np.array([[0, p[2], -p[1]],
                       [-p[2], 0, p[0]],
                       [p[1], p[0], 0]])

    return matrix

