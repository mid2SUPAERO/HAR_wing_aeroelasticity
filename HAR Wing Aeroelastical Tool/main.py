import numpy as np
from math import sin, cos, pi
from scipy.linalg import block_diag
from scipy.optimize import fsolve
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA

from plotaest3d import plotaest3d
from dinamicaflex import dinamicaflex
from Airplane import Airplane
from Structure import Structure
from Fuselage import Fuselage
from Engine import Engine
from Element import Element
from Node import Node
from flutter_speed import flutter_speed
from changedatarate import changedatarate


class Rot:
    def __init__(self, dihedral, sweep, twist):
        self.dihedral = dihedral
        self.sweep = sweep
        self.twist = twist


class AeroParams:
    def __init__(self, b, N, a, alpha0, clalpha, cm0, cd0, ndelta, cldelta, cmdelta):
        self.b = b
        self.N = N
        self.a = a
        self.alpha0 = alpha0
        self.clalpha = clalpha
        self.cm0 = cm0
        self.cd0 = cd0
        self.ndelta = ndelta
        self.cldelta = cldelta
        self.cmdelta = cmdelta


def create_uniform_structure(pos_cg, rot0, length, Inertia, mcs, KG, CG, aero, geometry, num_elements):
    ds = length / num_elements
    rot = rot0

    structure = Structure(num_elements)

    for i in range(0, num_elements):
        node1 = Node(mcs, pos_cg, Inertia, ((i - 1) * ds) / 20, aero, geometry)
        node2 = Node(mcs, pos_cg, Inertia, (ds / 2 + (i - 1) * ds) / 20, aero, geometry)
        node3 = Node(mcs, pos_cg, Inertia, i * ds / 20, aero, geometry)

        elem = Element(node1, node2, node3, rot, ds, KG, CG)

        structure.addElem(elem)

        rot.dihedral = 0
        rot.sweep = 0
        rot.twist = 0

    return structure


def create_flexible_member(num_elements, damp_ratio):
    Length = 16

    K11 = 1e10
    K22 = 1e4
    K33 = 2e4
    K44 = 4e6
    KG = np.diag([K11, K22, K33, K44])

    CG = damp_ratio * np.diag([K11, K22, K33, K44])

    pos_cg = np.array([0, 0.3, 0])



    I22 = 0.0
    I33 = 0.1
    I11 = 0.1
    mcs = 0.75
    Inertia = np.diag([I11, I22, I33])

    # aerodynamic data
    c = 1
    b = c / 2
    N = 0
    a = 0.0
    alpha0 = -5 * np.pi / 180
    clalpha = 2 * np.pi
    cm0 = 0
    cd0 = 0.02
    ndelta = 1
    cldelta = 0.01
    cmdelta = -0.1

    geometry = np.array([a, b])

    aero_right = AeroParams(b, N, a, alpha0, clalpha, cm0, cd0, ndelta, cldelta, cmdelta)

    rot0_right = Rot(dihedral=0, sweep=0*pi/180, twist=1*pi/180)

    right_wing = create_uniform_structure(pos_cg, rot0_right, Length, Inertia, mcs, KG, CG, aero_right, geometry,
                                          num_elements)

    # aerodynamic data
    c = 1
    b = c / 2
    N = 0
    a = 0.0
    alpha0 = 5 * np.pi / 180
    clalpha = 2 * np.pi
    cm0 = 0
    cd0 = 0.02
    ndelta = 1
    cldelta = -0.01
    cmdelta = 0.1

    geometry = np.array([a, b])

    aero_left = AeroParams(b, N, a, alpha0, clalpha, cm0, cd0, ndelta, cldelta, cmdelta)

    rot0_left = Rot(dihedral=pi, sweep=0*pi/180, twist=-1*pi/180)

    left_wing = create_uniform_structure(pos_cg, rot0_left, Length, Inertia, mcs, KG, CG, aero_left, geometry,
                                         num_elements)

    return right_wing,left_wing


def load_structure(num_elem, damp_ratio):
    right_wing, left_wing = create_flexible_member(num_elem, damp_ratio)

    right_wing.elementsVector[0].setH0(np.array([[0], [-0.3], [0], [1], [0], [0], [0], [1], [0], [0], [0], [1]]))
    left_wing.elementsVector[0].setH0(np.array([[0], [-0.3], [0], [1], [0], [0], [0], [1], [0], [0], [0], [1]]))
    right_wing.update()
    left_wing.update()

    fus = Fuselage(m=10, pcm=np.array([0, 0, 0]), I=np.zeros((3, 3)) + 0 * np.array([[0.2**2*10, 0, 0], [0, 0, 0], [0, 0, 0.2**2*10]]))

    Fmax = 1
    V0 = 1
    rho0 = 1
    nv = -1
    nrho = 0
    alphaf = 0
    betaf = 0
    numPI = 1
    numMEMB = 1
    numELEM = 1
    numNODE = 1

    motor1 = Engine(numPI, numMEMB, numELEM, numNODE, Fmax, rho0, V0, nrho, nv, alphaf, betaf)

    airp = Airplane(np.array([right_wing, left_wing]), np.array([motor1]), fus)

    # flexible_member = create_flexible_member(numele, damp_ratio)
    #
    # flexible_member.elementsVector[0].setH0(np.array([[0], [-0.0], [0], [1], [0], [0], [0], [1], [0], [0], [0], [1]]))
    # flexible_member.update()
    #
    # fus = Fuselage(m=0, pcm=np.array([0, 0, 0]), I=np.zeros((3, 3)) + 0 * np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])) # no fuselage
    # motor1 = Engine(0, 0, 0, 0) # no engines
    # airp = Airplane(np.array([flexible_member]), np.array([motor1]), fus)

    return airp


if __name__ == '__main__':
    # Airplane Initialization

    numele = 3
    damping = 0.04
    ap = load_structure(numele, damping)

    # Finds Equilibrium Condition
    altitude = 19931.7
    V = 10
    throttle = 0
    deltaflap = 0
    Vwind = 0

    isPinned = 'False'
    freeDEG = np.array([[1, 1, 1, 1, 1, 1]])

    rb_eq, strain_eq = ap.trimairplane(V, altitude, Vwind, throttle, deltaflap,freeDEG, isPinned)

    theta = rb_eq[0]
    deltaflap = rb_eq[1]
    engine_position = rb_eq[2]

    manete = rb_eq[2]
    thetaeq = theta
    straineq = strain_eq
    betaeq = np.array([[0, V * cos(thetaeq), - V * sin(thetaeq), 0, 0, 0]]).T
    keq = np.array([[thetaeq, 0, 0, altitude]]).T

    # ap.update(strain_eq * 0, np.zeros(strain_eq.shape), np.zeros(strain_eq.shape), np.zeros((int(np.sum(ap.membNAEDtotal)), 1)))
    # ap.plotAirplane3D()
    # tip_displacement = ap.members[0].elementsVector[numele-1].node3.h[2]
    #
    # print(tip_displacement)

    ap.update(strain_eq, np.zeros(strain_eq.shape), np.zeros(strain_eq.shape), np.zeros((int(np.sum(ap.membNAEDtotal)), 1)))
    ap.plotAirplane3D()
    tip_displacement = ap.members[0].elementsVector[numele-1].node3.h[2]

    print(tip_displacement)

# flutter speed - undeformed(Linear)

    # throttle = 0
    # deltaflap = 0
    # betaeq = np.zeros((6, 1))
    # keq = np.array([[0], [0], [0], [altitude]])
#     flut_speed = flutter_speed(20, 35, 0.01, ap, strain_eq * 0, altitude, betaeq,keq,manete,deltaflap, freeDEG)
#
#     print('flutter speed:')
#     print(flut_speed)
#     #print('flutter eig val:')
#     #print(flut_eig_val)
#
    # flutter speed - deformed(Nonlinear)
    flut_speed = flutter_speed(20, 100, 0.01, ap, strain_eq, altitude, betaeq,keq,manete,deltaflap, freeDEG)

    print('flutter speed:')
    print(flut_speed)
    #print('flutter eig val:')
    #print(flut_eig_val)

    # # Nonlinear simulation
    # # Initial conditions
    #
    # tSim = 5  # s
    # T = np.array([0.00, 0.49, 0.50, 1.99, 2.00, 3.49, 3.5, 100])
    # elev = np.array([0, 0, 1, 1, -1, -1, 0, 0]) / 5
    #
    # beta0 = np.array([[0], [V * cos(theta)], [-V * sin(theta)], [0], [0], [0]])
    # k0 = np.array([[theta], [0], [0], [altitude]])
    # strain0 = strain_eq
    # Vwind = 0
    #
    # tNL, strainNL, straindNL, lambdaNL, betaNL, kineticNL = ap.simulate(np.array([0, tSim]), strain0, beta0, k0, Vwind,
    #                                                                     engine_position, deltaflap, T, elev, freeDEG)
    #
    # dt = 0.1
    # ts, Xs = changedatarate(tNL, strainNL, dt)
    # ts, kinetics = changedatarate(tNL, kineticNL, dt)
    # tip_displacement = np.zeros((ts.shape[0], 1))
    #
    # for i in range(0, ts.shape[0]):
    #     ap.update(Xs[i, :], np.zeros(Xs[i, :].shape), np.zeros(Xs[i, :].shape), np.zeros((np.sum(ap.membNAEDtotal), 1)))
    #     tip_displacement[i] = ap.members[i].elementsVector[numele].node3.h[2]
    #
    # fig, ax = plt.subplots()
    # ax.plot(ts, tip_displacement, 'r')
    #
    # ax.set(xlabel='Time (s)', ylabel='Tip displacement (m)',
    #        title='Wing tip displacement')
    # ax.grid()
    #
    # longFig, longAx = plt.subplots(nrows=2, ncols=2)
    #
    # longAx[0, 0].plot(tNL, betaNL[:, 1], 'r')
    # longAx[0, 0].set_xlabel('t')
    # longAx[0, 0].set_ylabel('v (m/s)')
    #
    # longAx[0, 1].plot(tNL, betaNL[:, 2], 'r')
    # longAx[0, 1].set_xlabel('t')
    # longAx[0, 1].set_ylabel('w (m/s)')
    #
    # longAx[1, 0].plot(tNL, betaNL[:, 3], 'r')
    # longAx[1, 0].set_xlabel('t')
    # longAx[1, 0].set_ylabel('q (m/s)')
    #
    # longAx[1, 1].plot(tNL, kineticNL[:, 3], 'r')
    # longAx[1, 1].set_xlabel('t')
    # longAx[1, 1].set_ylabel('Altitude (m)')

    plt.show()

    # TODO ap.airplanemovie( ts', Xs,kinetics,dt,'test','gif')
