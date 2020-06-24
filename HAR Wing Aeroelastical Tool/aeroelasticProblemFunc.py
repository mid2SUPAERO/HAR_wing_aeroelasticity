import numpy as np
from math import sin, cos, pi
import matplotlib.pyplot as plt

from Airplane import Airplane
from Structure import Structure
from Fuselage import Fuselage
from Engine import Engine
from Element import Element
from Node import Node
from flutter_speed import flutter_speed


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


def create_flexible_member(num_elements, damp_ratio, rho, E, G, c, t, L, rootSweepAng, rootTwistAng):
    Length = L

    A = c * t  # Area
    Izz = t * c ** 3 / 12  # Area inertia z
    Iyy = c * t ** 3 / 12  # Area inertia y
    J = 0.33 * t * c ** 3

    K11 = E * A
    K22 = G * J
    K33 = E * Iyy
    K44 = E * Izz
    KG = np.diag([K11, K22, K33, K44])

    CG = damp_ratio * np.diag([K11, K22, K33, K44])

    pos_cg = np.array([0, 0, 0])

    I33 = rho * A * Izz
    I22 = rho * A * Iyy
    I11 = I22 + I33
    mcs = rho * A
    Inertia = np.diag([I11, I22, I33])

    # aerodynamic data
    b = c / 2
    N = 4
    a = 0
    cm0 = 0
    cd0 = 0.02
    ndelta = 0

    # for the right wing
    alpha0 = -5 * np.pi / 180
    clalpha = 2 * np.pi
    cldelta = 0
    cmdelta = 0

    geometry = np.array([a, b])

    aero_right = AeroParams(b, N, a, alpha0, clalpha, cm0, cd0, ndelta, cldelta, cmdelta)

    rot0_right = Rot(dihedral=0, sweep=rootSweepAng, twist=rootTwistAng)

    right_wing = create_uniform_structure(pos_cg, rot0_right, Length, Inertia, mcs, KG, CG, aero_right, geometry,
                                          num_elements)

    return right_wing


def load_structure(num_elem, damp_ratio, rho, E, G, c, t, L, rootSweepAng, rootTwistAng):
    flexible_member = create_flexible_member(num_elem, damp_ratio, rho, E, G, c, t, L, rootSweepAng, rootTwistAng)

    flexible_member.elementsVector[0].setH0(np.array([[0], [-0.0], [0], [1], [0], [0], [0], [1], [0], [0], [0], [1]]))
    flexible_member.update()

    fus = Fuselage(m=0, pcm=np.array([0, 0, 0]),
                   I=np.zeros((3, 3)) + 0 * np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]))  # no fuselage
    motor1 = Engine(0, 0, 0, 0)  # no engines
    airp = Airplane(np.array([flexible_member]), np.array([motor1]), fus)

    return airp


def aeroelasticProblemFunc(t, c, L, rootSweepAng, rootTwistAng):
    # ---------- Material Properties----------#
    # AL 7075
    rho = 2810  # kg/m³
    E = 71.7e9  # Pa
    G = 29.9e9  # Pa

    mass = rho * c[0] * t[0] * L[0]

    # Airplane Initialization
    numele = 8
    damping = 0.04
    ap = load_structure(numele, damping, rho, E, G, c[0], t[0], L[0], rootSweepAng[0], rootTwistAng[0])

    # Finds Equilibrium Condition
    altitude = 19931.7
    V = 0
    throttle = 0
    deltaflap = 0
    Vwind = 10

    isPinned = 'True'
    freeDEG = np.array([[0, 0, 0, 0, 0, 0]])

    rb_eq, strain_eq = ap.trimairplane(V, altitude, Vwind, throttle, deltaflap, freeDEG, isPinned)

    theta = rb_eq[0]
    deltaflap = rb_eq[1]
    engine_position = rb_eq[2]

    manete = rb_eq[2]
    thetaeq = theta
    straineq = strain_eq
    betaeq = np.array([[0, V * cos(thetaeq), - V * sin(thetaeq), 0, 0, 0]]).T
    keq = np.array([[thetaeq, 0, 0, altitude]]).T

    ap.update(strain_eq, np.zeros(strain_eq.shape), np.zeros(strain_eq.shape),
              np.zeros((int(np.sum(ap.membNAEDtotal)), 1)))

    ap.plotAirplane3D()

    tip_displacement = ap.members[0].elementsVector[numele - 1].node3.h[2]

    # flutter speed - deformed(Nonlinear)
    flut_speed = flutter_speed(10, 100, 0.1, ap, strain_eq, altitude, betaeq, keq, manete, deltaflap, freeDEG)

    J = -flut_speed / 10 + mass / 35.125

    plt.show()

    return J, flut_speed, mass, tip_displacement[0]
