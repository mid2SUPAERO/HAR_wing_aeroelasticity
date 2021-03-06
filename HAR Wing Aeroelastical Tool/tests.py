from AeroParamsUnsteady import AeroParamsUnsteady
from Node import Node
import math
import numpy as np
import matplotlib.pyplot as plt
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

'''
aeroParams = AeroParamsUnsteady(2, 0.5, 0, -5 * math.pi / 180, 2 * math.pi, 0.01, 0, -0.1, 0.02, 1)

print('E1 = ')
print(aeroParams.E1)
print('E2 = ')
print(aeroParams.E2)
print('E3 = ')
print(aeroParams.E3)
print('E4 = ')
print(aeroParams.E4)
print(aeroParams.A)
print(aeroParams.B)
print(aeroParams.clAlpha)

print('-----noh------')
inertia = np.array([[0,0,0],[0,0.1,0],[0,0,0.1]])
cg = np.array([0,0.3,0])

noh = Node(0.75, cg, inertia, 0.2, aeroParams)
print(noh.aeroParams.clAlpha)
print(noh.massMatrix)

noh.aeroParams = aero = AeroParamsUnsteady(2, 0.5, 0, -5 * math.pi / 180, 2 * math.pi, 0.01, 0, -0.1, 0.02, 1)
print(aero.B)
print(noh.aeroParams.B)

aero.B = np.array([[3], [4]])

print(aero.B)
print(aeroParams.B)
print(noh.aeroParams.B)
print('-----element------')
'''

#
# def inp(x):
#     return x


'''def residual(t, y, yd):
    res_0 = yd[0] - y[2]
    res_1 = yd[1] - y[3]
    res_2 = yd[2] + y[4] * y[0]
    res_3 = yd[3] + y[4] * y[1] + 9.82
    res_4 = y[2] ** 2 + y[3] ** 2 - y[4] * (y[0] ** 2 + y[1] ** 2) - y[1] * 9.82

    return np.array([res_0, res_1, res_2, res_3, res_4])


t0 = 0.0  # Initial time
y0 = [1.0, 0.0, 0.0, 0.0, 0.0]  # Initial states
yd0 = [0.0, 0.0, 0.0, -9.82, 0.0]  # Initial state derivatives

model = Implicit_Problem(residual, y0, yd0, t0)  # Create an Assimulo problem
model.name = 'Pendulum'  # Specifies the name of problem (optional)

sim = IDA(model)
tfinal = 10.0  # Specify the final time
ncp = 500  # Number of communication points (number of return points)

t, y, yd = sim.simulate(tfinal,
                        ncp)  # Use the .simulate method to simulate and provide the final time and ncp (optional)
# Plot the result
sim.plot()'''

# def weissinger(t, y, yp, test1):
#     return t * y ** test1(2) * yp ** 3 - y ** 3 * yp ** 2 + t * (t ** 2 + 1) * yp - t ** 2 * y
#
#
# t0 = 1
# y0 = math.sqrt(3 / 2)
# yp0 = 0.8165
#
# model2 = Implicit_Problem(weissinger, y0, yp0, t0, p0=inp())
# model2.name = 'Weissinger'
# sim = IDA(model2)
#
# tfinal = 10
# ncp = 100
#
# t, y, yp = sim.simulate(tfinal, ncp)
#
# sim.plot()

'''t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2 * np.pi * t)


fig, ax = plt.subplots()
ax.plot(t, s,'r')

ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       title='About as simple as it gets, folks')
ax.grid()

fig2, ax2 = plt.subplots()
ax2.plot(s, t,'g')

ax2.set(xlabel='time (s)', ylabel='voltage (mV)',
        title='About as simple as it gets, folks')
ax2.grid()

plt.show()'''

'''def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

ax = plt.axes(projection='3d')

ax.set_title('surface')

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='viridis')

plt.show()'''
from aeroelasticProblemFunc import aeroelasticProblemFunc

t = np.array([0.0010001])
c = np.array([0.3911399])
L = np.array([4.80474242])
sweep = np.array([0.04876023])
twist = np.array([0.01767856])

J, flut_speed, mass, tip_displacement = aeroelasticProblemFunc(t, c, L, sweep, twist)

print(J)
print(flut_speed)
print(mass)
print(tip_displacement)
