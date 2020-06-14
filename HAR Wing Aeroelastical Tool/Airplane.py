import numpy as np
import math
from scipy.linalg import block_diag
from scipy.optimize import fsolve
from plotaest3d import plotaest3d
from dinamicafex import dinamicaflex
from Fuselage import Fuselage


def veclin(i, numstates):
    vec = np.zeros((1, numstates))
    vec[0, i] = 1

    return vec


def equilibracorpo(vec, strain, strainp, strainpp, lambd, beta, betap, kinetic, ap, V):
    theta = vec[0]
    deltaflap = vec[1]
    tracao = vec[2]
    beta[2] = -V * math.sin(theta)
    beta[1] = V * math.cos(theta)
    Vento = 0
    kinetic[1] = theta

    FLAG = 0
    Xp, bp, lambdap, kp = dinamicaflex(0, strain, strainp, strainpp, lambd, beta, betap, kinetic, ap, Vento, tracao,
                                       np.array([deltaflap, 0, 0]), FLAG)
    zero = 10000 * np.array([bp[2], bp[3], bp[4]]).transpose()

    return zero


def equilibraestrutura(strain, strainp, strainpp, lambd, beta, betap, kinetic, ap, V, tracao, deltaflap):
    FLAG = 0
    Xp, bp, lambdap = dinamicaflex(0, strain, strainp, strainpp,
                                   lambd, beta, betap, kinetic, ap, V, tracao, deltaflap, FLAG)
    zero = Xp

    return zero


def equilibratudo(vec, strain, strainp, strainpp, lambd, beta, betap, kinetic, ap, V, tracao, deltaflap):
    # global softPARAMS # TODO necessario?
    strain = vec[0:strainp.shape[1]]
    theta = vec[strainp.shape[1] + 1]
    deltaflap = vec[strainp.shape[1] + 2]
    tracao = vec[strainp.shape[1] + 3]  # TODO por que?
    beta[3] = -V * math.sin(theta)
    beta[2] = V * math.cos(theta)
    Vento = 0
    kinetic[1] = theta

    # softPARAMS.isEQ = 1 # TODO masoq???
    FLAG = 0
    Xp, bp, lambdap = dinamicaflex(0, strain, strainp, strainpp,
                                   lambd, beta, betap, kinetic, ap, V, tracao, deltaflap, FLAG)
    # softPARAMS.isEQ = 0
    zero = np.array([[Xp], [bp[2]], [bp[3]], [bp[4]]])

    return zero


class Airplane:
    def __init__(self, structures, engines, fus=Fuselage(0, np.array([0, 0, 0]), np.zeros((3, 3)))):
        self.NUMmembers = structures.shape[1]
        self.members = structures
        self.NUMele = 0
        self.prop = engines
        self.membSIZES = np.zeros((1, self.NUMmembers))
        self.membNAED = np.zeros((1, self.NUMmembers))
        self.membNAEDtotal = np.zeros((1, self.NUMmembers))
        self.fus = fus

        self.Jhep = None
        self.Jpep = None
        self.Jthetaep = None
        self.Jhb = None
        self.Jpb = None
        self.Jthetab = None

        for ii in range(0, self.NUMmembers):
            self.membSIZES[ii] = self.members[ii].numOfElems
            self.membNAED[ii] = self.members[ii].elementsVector[1].node1.aeroParams.N

            self.membNAEDtotal[ii] = 3 * self.membNAED[ii] * self.membSIZES[ii]
            self.NUMele = self.NUMele + self.membSIZES[ii]

        self.NUMaedstates = np.sum(self.membNAEDtotal)

        for ii in range(0, self.prop.shape[1]):
            self.prop[ii].NODEpos = 1 + 3 * np.sum(self.membSIZES[0:(self.prop[ii].numMEMB - 1)]) + (
                    self.prop[ii].numELEM - 1) * 3 + (self.prop[ii].numNODE - 1)

        Me = np.array([])
        K = np.array([])
        C = np.array([])
        N = np.array([])
        B = np.array([])

        for ii in range(0, self.NUMmembers):
            Me = block_diag(Me, self.members[ii].getMe())
            K = block_diag(K, self.members[ii].getK())
            C = block_diag(C, self.members[ii].getC())
            N = np.block([[N], [self.members[ii].getN()]])
            B = block_diag(B, self.members[ii].getB())

        self.Me = Me
        self.K = K
        self.C = C
        self.N = N
        self.B = B

    # ------------- # ------------- METHODS ------------- # ------------- #

    def updateStrJac(self):
        Jhep = []
        Jpep = []
        Jthetaep = []
        Jhb = []
        Jpb = []
        Jthetab = []

        for i in range(0, self.NUMmembers):
            self.members[i].elementsVector[1].memberJhep = self.members[i].getJhep()
            Jhep = block_diag(Jhep, self.members[i].elementsVector[1].memberJhep)

            self.members[i].elementsVector[1].memberJpep = self.members[i].getJpep(
                self.members[i].elementsVector[1].memberJhep)
            Jpep = block_diag(Jpep, self.members[i].elementsVector[1].memberJpep)

            self.members[i].elementsVector[1].memberJthetaep = self.members[i].getJthetaep(
                self.members[i].elementsVector, self.members[i].elementsVector[1].memberJhep)  # TODO checar
            Jthetaep = block_diag(Jthetaep, self.members[i].elementsVector[1].memberJthetaep)

            self.members[i].elementsVector[1].memberJhb = self.members[i].getJhb()
            Jhb = np.block([[Jhb], [self.members[i].elementsVector[1].memberJhb]])

            self.members[i].elementsVector[1].memberJpb = self.members[i].getJpb(
                self.members[i].elementsVector[1].memberJhb)
            Jpb = [[Jpb], [self.members[i].elementsVector[1].memberJpb]]

            self.members[i].elementsVector[1].memberJthetab = self.members[i].getJthetab()
            Jthetab = [[Jthetab], [self.members[i].elementsVector[1].memberJthetab]]

        self.Jhep = Jhep
        self.Jpep = Jpep
        # if isempty(self.Jthetaep) # TODO talvez tenha que adicionar depois? == None ou == []
        self.Jthetaep = Jthetaep
        self.Jhb = Jhb
        self.Jpb = Jpb
        self.Jthetab = Jthetab

    def plotAirplane3D(self, translate=np.array([0, 0, 0])):
        # hold on # TODO ajustar para plot no python
        ploti = np.zeros((1, self.NUMmembers))

        for i in range(0, self.NUMmembers):
            ploti[i] = plotaest3d(self.members[i], translate)

        return ploti

    def update(self, strain, strainp, strainpp, lambd):
        for i in range(0, self.NUMmembers):
            self.members[i].elementsVector[1].strainm = strain[(np.sum(self.membSIZES[0:(i - 1)]) * 4):(
                    np.sum(self.membSIZES[0:i]) * 4)].tanspose()
            self.members[i].elementsVector[1].strainpm = strainp[(np.sum(self.membSIZES[0:(i - 1)]) * 4):(
                    np.sum(self.membSIZES[0:i]) * 4)].tanspose()
            self.members[i].elementsVector[1].strainppm = strainpp[(np.sum(self.membSIZES[0:(i - 1)]) * 4):(
                    np.sum(self.membSIZES[0:i]) * 4)].tanspose()
            self.members[i].elementsVector[1].lambdm = lambd[(np.sum(self.membNAEDtotal[0:(i - 1)])):(
                np.sum(self.membNAEDtotal[0:i]))].tanspose()

        for i in range(0, self.NUMmembers):
            for j in range(0, self.membSIZES[i]):
                self.members[i].elementsVector[j].setStrain(
                    self.members[i].elementsVector[1].strainm[(j * 4): (4 + j * 4)],
                    self.members[i].elementsVector[1].strainpm[(j * 4): (4 + j * 4)])

        for i in range(0, self.NUMmembers):
            self.members[i].update()

    def airplanemovie(self):  # ah mas vai pass ar ainda por muito tempo
        pass

    def dinamicaflexODEimplicit(self): # TODO resolver
        pass

    def dinamicaflexODE(self):
        pass


    def linearize(self, straineq, betaeq, keq, manete, deltaflap, Vwind):

        if self.NUMaedstates == 0:
            lambda0 = []
        else:
            lambda0 = np.zeros((self.NUMaedstates, 1))

        Airplane.update(self, straineq, straineq * 0, straineq * 0, lambda0)
        Airplane.updateStrJac(self)

        x = np.block([[straineq.tanspose()], [np.zeros(np.shape(straineq.transpose()))], [lambda0], [betaeq], [keq]])
        xp = np.zeros((np.shape(x), np.shape(x)))

        delta = 1e-9

        Alin = np.zeros(
            ((self.NUMaedstates + straineq.shape[1] * 2 + 10), (self.NUMaedstates + np.shape(straineq, 2) * 2 + 10)))
        Mlin = np.eye((self.NUMaedstates + straineq.shape[1] * 2 + 10))

        for i in range(0, xp.shape[0]):
            soma = Airplane.dinamicaflexODEimplicit(self, 0, x + veclin(i, x.shape[0]).tanspose() * delta, xp, Vwind,
                                                    manete, deltaflap)
            subtr = Airplane.dinamicaflexODEimplicit(self, 0, x - veclin(i, x.shape[0]).tanspose() * delta, xp, Vwind,
                                                     manete, deltaflap)

            if (i >= (straineq.shape[1] + 1)) and (i <= (straineq.shape[1] * 2)) or (
                    (i >= self.NUMaedstates + straineq.shape[1] * 2 + 1) and (
                    i <= self.NUMaedstates + straineq.shape[1] * 2 + 6)):  # TODO talvez precise alterar esse i
                somap = Airplane.dinamicaflexODEimplicit(self, 0, x, xp + veclin(i, x.shape[0]).tanspose() * delta,
                                                         Vwind, manete, deltaflap)
                subtrp = Airplane.dinamicaflexODEimplicit(self, 0, x, xp - veclin(i, x.shape[0]).tanspose() * delta,
                                                          Vwind, manete, deltaflap)
                Mlin[:, i] = -(somap - subtrp) / (2 * delta)

            Alin[:, i] = (soma - subtr) / (2 * delta)

        A = np.linalg.solve(Mlin, Alin)

        Aaero = np.linalg.solve(
            Mlin[0:(self.NUMaedstates + straineq.shape[1] * 2), 0: (self.NUMaedstates + straineq.shape[1] * 2)],
            Alin[0: (self.NUMaedstates + straineq.shape[1] * 2), 0: (self.NUMaedstates + straineq.shape[1] * 2)])
        Abody = np.linalg.solve(
            Mlin[(self.NUMaedstates + straineq.shape[1] * 2):(self.NUMaedstates + straineq.shape[1] * 2 + 10),
            (self.NUMaedstates + straineq.shape[1] * 2): (self.NUMaedstates + straineq.shape[1] * 2 + 10)],
            Alin[(self.NUMaedstates + straineq.shape[1] * 2): (self.NUMaedstates + straineq.shape[1] * 2 + 10),
            (self.NUMaedstates + straineq.shape[1] * 2): (self.NUMaedstates + straineq.shape[1] * 2 + 10)])

        return A, Aaero, Abody

    def trimairplane(self, V, H, Vwind, tracao=0, deltaflap=0):
        strain = np.zeros((self.NUMele * 4, 1)).transpose()
        strainp = np.zeros((self.NUMele * 4, 1)).transpose()
        strainpp = np.zeros((self.NUMele * 4, 1)).transpose()

        lambd = np.zeros((np.sum(self.membNAEDtotal), 1))
        Airplane.updateStrJac(self)

        beta = np.zeros((6, 1))
        betap = np.zeros((6, 1))

        theta = 0
        phi = 0
        psi = 0
        kinetic = np.array([theta, phi, psi, H])

        # TODO não entendo por que ele resolve duas vezes, talvez pra plotar, mas deve dar pra fazer 1 vez so, verificar dps
        # TODO sinto que vai dar problema aqui
        vec = fsolve(equilibracorpo, np.array([0, 0, 0]),
                     args=(strain, strainp, strainpp, lambd, beta, betap, kinetic, self, V,),
                     xtol=1e-15, maxfev=20000)
        theta = vec[0]
        kinetic[0] = theta
        deltaflap = vec[1]
        tracao = vec[2]
        beta[2] = -V * math.sin(theta)
        beta[1] = V * math.cos(theta)
        strainEQ = fsolve(equilibraestrutura, strain,
                          args=(strainp, strainpp, lambd, beta, betap, kinetic, self, Vwind, tracao,
                                np.array([deltaflap, 0, 0])),
                          xtol=1e-15, maxfev=20000)

        # figure(100); # TODO adaptar
        # plotairplane3d(ap);

        theta = 0
        vec = np.array([0, 0, 0])
        kinetic[0] = theta
        beta[2] = -0 * math.sin(theta)
        beta[1] = 0 * math.cos(theta)
        strainEQ = fsolve(equilibraestrutura, strain,
                          args=(strainp, strainpp, lambd, beta, betap, kinetic, self, Vwind, tracao, deltaflap),
                          xtol=1e-17, maxfev=20000)

        # figure(100); # TODO adaptar
        # plotairplane3d(ap);

        # TODO aqui ele printava esses valores, vale a pena? precisa? se não, tirar
        # theta
        # deltaflap
        # tracao

        return vec, strainEQ

    def trimairplanefull(self):  # TODO necessario? prefiro usar so esse ou o full
        pass

    def simulate(self, tspan, strain0, beta0, kinetic0, Vwind, manete, deltaflap):
        x0 = np.block(
            [[strain0.transpose()], [np.zeros((self.NUMele * 4, 1))], [np.zeros((self.NUMaedstates, 1))], [beta0],
             [kinetic0]])

        xp0 = np.zeros(x0.shape)

        # options = odeset('OutputFcn', @ odeprog, 'Events', 'MaxStep', 0.10)
        # TODO resolver a treta da ode
        # [t, X] = ode15i( @ (t, x, xp) dinamicaflexODEimplicit(t, x, xp, ap, Vwind, manete(t), deltaflap(t)), tspan, x0, xp0, options);

        strain = X[:, 0: self.NUMele * 4]
        strainp = X[:, self.NUMele * 4: 2 * self.NUMele * 4]
        lambd = X[:, 2 * self.NUMele * 4: (2 * self.NUMele * 4 + self.NUMaedstates)]
        beta = X[:, 2 * self.NUMele * 4 + self.NUMaedstates: 2 * self.NUMele * 4 + self.NUMaedstates + 6]
        kinetic = X[:, 2 * self.NUMele * 4 + self.NUMaedstates + 6: 2 * self.NUMele * 4 + self.NUMaedstates + 10]

        return t, strain, strainp, lambd, beta, kinetic

