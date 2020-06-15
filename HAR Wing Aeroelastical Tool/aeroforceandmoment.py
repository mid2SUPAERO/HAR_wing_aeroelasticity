import numpy as np
from math import sin, cos, atan, pi


def aeroforceandmoment(strain, strainp, strainpp, lambd, beta, betap, membro, Vwindy, rho, deltaflap):
    Vwind = np.array([[0], [Vwindy], [0]])
    membersize = membro.numOfElems
    pp = membro.elementsVector[1].memberJpep * strainp + membro.elementsVector[1].memberJpb * beta
    ppp = membro.elementsVector[1].memberJpep * strainpp + membro.elementsVector[1].memberJpb * betap
    thetap = membro.elementsVector[1].memberJthetaep * strainp + membro.elementsVector[1].memberJthetab * beta
    thetapp = membro.elementsVector[1].memberJthetaep * strainpp + membro.elementsVector[1].memberJthetab * betap
    tao = 15

    # TODO talvez de errado
    F = np.zeros((9 * (membersize + 1), 1))
    M = np.zeros((9 * (membersize + 1), 1))
    lambdap = np.zeros((((membersize + 1) * 3 * membro.elementsVector[1].node1.aeroParams.N), 1))

    for i in range(0, membersize):
        # Node 1
        N = membro.elementsVector[i].node1.aeroParams.N
        lambda0 = 0.5 * lambd[(i * N * 3):(N + i * N * 3)] @ membro.elementsVector[i].node1.aeroParams.B
        b = membro.elementsVector[i].node1.aeroParams.b
        d = membro.elementsVector[i].node1.aeroParams.a * b

        Forca, Momento, flambda = getforceandmoment(membro.elementsVector[i].node1.h, pp[(i * 9):(3 + i * 9)],
                                                    thetap[(i * 9): (3 + i * 9)], ppp[(i * 9): (3 + i * 9)],
                                                    thetapp[(i * 9): (3 + i * 9)], lambda0, Vwind, b, rho, d,
                                                    membro.elementsVector[i].node1.aeroParams,
                                                    lambd[(i * N * 3):(N + i * N * 3)], deltaflap)
        F[(i * 9): (3 + i * 9), 0] = Forca
        M[(i * 9): (3 + i * 9), 0] = Momento
        lambdap[(i * 3 * N): (N + i * 3 * N), 0] = flambda

        # Node 2
        N = membro.elementsVector[i].node2.aeroParams.N
        lambda0 = 0.5 * lambd[(N + i * N * 3):(2 * N + i * N * 3)] @ membro.elementsVector[i].node2.aeroParams.B
        b = membro.elementsVector[i].node2.aeroParams.b
        d = membro.elementsVector[i].node2.aeroParams.a * b

        Forca, Momento, flambda = getforceandmoment(membro.elementsVector[i].node2.h, pp[(3 + i * 9):(6 + i * 9)],
                                                    thetap[(3 + i * 9): (6 + i * 9)], ppp[(3 + i * 9): (6 + i * 9)],
                                                    thetapp[(3 + i * 9): (6 + i * 9)], lambda0, Vwind, b, rho, d,
                                                    membro.elementsVector[i].node2.aeroParams,
                                                    lambd[(N + i * N * 3):(2 * N + i * N * 3)], deltaflap)
        F[(3 + i * 9): (6 + i * 9), 0] = Forca
        M[(3 + i * 9): (6 + i * 9), 0] = Momento
        lambdap[(N + i * 3 * N): (2 * N + i * 3 * N), 0] = flambda

        # Node 3
        N = membro.elementsVector[i].node3.aeroParams.N
        lambda0 = 0.5 * lambd[(2 * N + i * N * 3):(3 * N + i * N * 3)] @ membro.elementsVector[i].node3.aeroParams.B
        b = membro.elementsVector[i].node3.aeroParams.b
        d = membro.elementsVector[i].node3.aeroParams.a * b

        Forca, Momento, flambda = getforceandmoment(membro.elementsVector[i].node3.h, pp[(6 + i * 9):(9 + i * 9)],
                                                    thetap[(6 + i * 9): (9 + i * 9)], ppp[(6 + i * 9): (9 + i * 9)],
                                                    thetapp[(6 + i * 9): (9 + i * 9)], lambda0, Vwind, b, rho, d,
                                                    membro.elementsVector[i].node3.aeroParams,
                                                    lambd[(2 * N + i * N * 3):(3 * N + i * N * 3)], deltaflap)
        F[(6 + i * 9): (9 + i * 9), 0] = Forca
        M[(6 + i * 9): (9 + i * 9), 0] = Momento
        lambdap[(2 * N + i * 3 * N): (3 * N + i * 3 * N), 0] = flambda

    return F, M, lambdap


def getforceandmoment(h, dotp, dottheta, ddotp, ddottheta, lambda0, Vwind, b, rho, d, aero, lambd, deltaflap):
    e1 = np.array([1, 0, 0]).transpose()
    e2 = np.array([0, 1, 0]).transpose()
    e3 = np.array([0, 0, 1]).transpose()

    if aero.ndelta > 0:
        delta = deltaflap[aero.ndelta - 1]  # TODO verificar o -1
    else:
        delta = 0

    CBw = np.block([h[3:6], h[6: 9], h[9: 12]])
    Cwa0 = np.eye(3)
    alpha0 = aero.alpha0
    Cwa0 = np.block([[1, 0, 0], [0, cos(alpha0), sin(alpha0)], [0, - sin(alpha0), cos(alpha0)]])

    CBa0 = CBw @ Cwa0
    doty = e2.transpose() @ CBa0.transpose() @ (dotp + Vwind)  # TODO talvez não seja @ no segundo
    dotz = e3.transpose() @ CBa0.transpose() @ (dotp + Vwind)
    alphat = atan(-dotz / doty)
    Ca0a1 = np.block([[1, 0, 0], [0, cos(alphat), sin(alphat)], [0, - sin(alphat), cos(alphat)]])

    CBA = CBa0 @ Ca0a1

    ddotz = e3.transpose() @ CBa0.transpose() @ ddotp
    dotalpha = e1.transpose() @ CBa0.transpose() @ dottheta
    ddotalpha = e1.transpose() @ CBa0.transpose() @ ddottheta

    clalpha = aero.clalpha
    cldelta = aero.cldelta
    cd0 = aero.cd0
    cm0 = aero.cm0
    cmdelta = aero.cmdelta

    # TODO todas as variaveis abaixo xão escalares?
    L = pi * rho * b ^ 2 * (-ddotz + doty * dotalpha - d * ddotalpha) + clalpha * rho * b * doty ^ 2 * (
                -dotz / doty + (b / 2 - d) * dotalpha / doty - lambda0 / doty) + rho * b * doty ^ 2 * cldelta * delta
    M = -pi * rho * b ^ 3 * (
                -0.5 * ddotz + doty * dotalpha + (0.125 * b - 0.5 * d) * ddotalpha) + 2 * rho * b ^ 2 * doty ^ 2 * (
                    cm0 + cmdelta * delta)
    D = -rho * b * doty ^ 2 * cd0

    Forca = CBA @ np.array([[0], [D], [L]])
    Mx = M + (d + 0.5 * b) * (L * cos(alphat + alpha0) - D * sin(alphat + alpha0))
    Momento = CBA * np.array([[Mx], [0], [0]])

    Vwind = doty
    # TODO checar se tem mais @
    flambda = aero.E1 @ lambd.transpose() * np.linalg.norm(Vwind) / b + aero.E2 * (
                -ddotz / b) + aero.E3 * ddotalpha + aero.E4 * dotalpha * np.linalg.norm(Vwind) / b

    return Forca, Momento, flambda
