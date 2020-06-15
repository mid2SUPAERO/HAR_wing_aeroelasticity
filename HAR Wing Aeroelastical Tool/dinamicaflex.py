import numpy as np
import math

from matrixcross import matrixcross
from OMEGA import OMEGA
from atmosfera import atmosfera
from getJhbp import getJhpb
from aeroforceandmoment import aeroforceandmoment


# TODO matlab precisa colocar o tempo, mas aqui é necessario?
def dinamicaflex(t, strain, strainp, strainpp, lambd, beta, betap, kinetic, ap, V, manete, deltaflap, FLAG):
    H = kinetic[4]
    isEQ = None
    isIMPLICIT = None

    g = 9.80665  # m/s**2 standard earth gravity
    rho = atmosfera(H)

    # TODO verificar se todos são necessarios
    if FLAG == 0:
        isEQ = 1

    elif FLAG == 1:
        isEQ = 0
        isIMPLICIT = 1

    elif FLAG == 2:
        isEQ = 0
        isIMPLICIT = 0

    ap.update(strain, strainp, strainpp, lambd)

    if isEQ == 1:
        ap.updateStrJac()

    KFF = ap.K

    if isEQ == 0:
        MFF = ap.Jhep.transpose() @ ap.Me @ ap.Jhep
        CFF = ap.C
        MFB = ap.Jhep.transpose() @ ap.Me @ ap.Jhb
        MBF = ap.Jhb.transpose() @ ap.Me @ ap.Jhep
        MRB = ap.fus.MRB
        MBB = ap.Jhb.transpose() @ ap.Me @ ap.Jhb + MRB

        CBF = ap.Jhb.transpose() @ ap.Me @ ap.Jhep * 0
        omega = matrixcross(beta[3:6]).transpose()
        nodesnumber = ap.NUMele * 3
        Hhb = OMEGA(omega, nodesnumber) @ ap.Jhb

        hp = ap.Jhep @ strainp.transpose()
        Jhbp = getJhpb(hp, nodesnumber)
        CFB = ap.Jhep.transpose() @ ap.Me @ Hhb + 2 * ap.Jhep.transpose() @ ap.Me @ Jhbp
        CRB = ap.fus.getCRB(omega)
        CBB = ap.Jhb.transpose() @ ap.Me @ Hhb + 2 * ap.Jhb.transpose() @ ap.Me @ Jhbp + CRB

    FAERO = np.array([])
    MAERO = np.array([])
    FLAMBDA = np.array([])

    for i in range(0, ap.NUMmembers):
        [FAEROm, MAEROm, FLAMBDAm] = aeroforceandmoment(ap.members[i].elementsVector[1].strainm,
                                                        ap.members[i].elementsVector[1].strainpm,
                                                        ap.members[i].elementsVector[1].strainppm,
                                                        ap.members[i].elementsVector[1].lambdam, beta, betap,
                                                        ap.members[i], V, rho, deltaflap)

        FAERO = np.block([[FAERO], [FAEROm]])
        MAERO = np.block([[MAERO], [MAEROm]])
        FLAMBDA = np.block([[FLAMBDA], [FLAMBDAm]])

    U = beta[2]

    FPROP = ap.prop.getFPROP(ap, manete, rho, U)

    theta = kinetic[1]
    phi = kinetic[2]
    psi = kinetic[3]
    H = kinetic[4]  # TODO de novo? necessario?

    GX = g * math.sin(phi) * math.cos(theta)
    GY = -g * math.sin(theta)
    GZ = -g * math.cos(phi) * math.cos(theta)
    GRAVITY = np.block([[GX], [GY], [GZ]])

    RAEROF = ap.Jpep.transpose() @ ap.B @ FAERO + ap.Jpep.transpose() @ FPROP + ap.Jthetae.transpose() @ ap.B @ MAERO + ap.Jhep.transpose() @ ap.N @ GRAVITY
    RAEROB = ap.Jpb.transpose() @ ap.B @ FAERO + ap.Jpb.transpose() @ FPROP + ap.Jthetab.transpose() @ ap.B @ MAERO + ap.Jhb.transpose() @ ap.N @ GRAVITY + ap.fus.N @ GRAVITY

    if isEQ == 0:
        if isIMPLICIT == 1:
            xp = np.linalg.solve(MFF, (
                    -MFB @ betap - CFB @ beta - CFF @ strainp.transpose() - KFF @ strain.transpose() + RAEROF))
            bp = np.linalg.solve(MBB, (
                    -MBF @ strainpp.transpose() - CBF @ strainp.transpose() - CBB @ beta + RAEROB) * np.array(
                [1, 1, 1, 1, 1, 1]).transpose())
        else:
            xpbp = np.linalg.solve(np.block([[MFF, MFB], [MBF, MBB]]), (
                    np.block([-CFF, -CFB], [-CBF, -CBB]) @ np.block([[strainp.transpose()], [beta]]) - np.block(
                [[KFF @ strain.transpose()], [np.zeros((6, 1))]]) + np.block([[RAEROF], [RAEROB]])))
            xp = xpbp[0:ap.NUMele * 4]
            bp = xpbp[(ap.NUMele * 4):(ap.NUMele * 4 + 6)] * np.linalg.solve.transpose()
    else:
        xp = - KFF @ strain.transpose() + RAEROF
        bp = + RAEROB * np.array([1, 1, 1, 1, 1, 1]).transpose()

    lambdap = FLAMBDA

    P = beta[4]
    Q = beta[3]
    R = -beta[5]
    thetap = Q * math.cos(phi) - R * math.sin(phi)
    phip = P + math.tan(theta) * (Q * math.sin(phi) + R * math.cos(phi))
    psip = (Q * math.sin(phi) + R * math.cos(phi)) / math.cos(theta)
    U = beta[1]
    V = beta[0]
    W = -beta[2]
    Hp = U * math.sin(theta) - V * math.sin(phi) * math.cos(theta) - W * math.cos(phi) * math.cos(theta)

    kineticp = np.array([[thetap], [phip], [psip], [Hp]])

    return xp, bp, lambdap, kineticp
