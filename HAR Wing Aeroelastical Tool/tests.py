from AeroParamsUnsteady import AeroParamsUnsteady
from Node import Node
import math
import numpy as np

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


print('-----element------')
