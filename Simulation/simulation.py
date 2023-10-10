'''This file houses the simulation.'''
import constants as c
from tools import kepler2rv
import cr3bp
from scipy import integrate as int
import numpy as np
import matplotlib.pyplot as plt


## Setup timestep of propagation
t0 = 0              # Initial tTime
tbound = 100000      # Final time
tsteps = 10000

## Create state vector
r, v = kepler2rv(200000*1000, .01, 30, 0, 100, 0, c.G*c.earthMass)
state = np.concatenate((r,v),axis=None)


## Pass in a state vector to a certain dynamics model   
# dynamics = cr3bp.SynodicEOMs(t, state, mu)

## Integrate those dynamics with a certian solver
odesol = int.solve_ivp(fun=cr3bp.SynodicEOMs, t_span=[t0,tbound], y0=state, method='DOP853', max_step=(tbound-t0)/tsteps)

# rvecs = []
# for state in odesol.y.T:

x_vals = odesol.y[0,:]
y_vals = odesol.y[1,:]
z_vals = odesol.y[2,:]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x_vals, y_vals, z_vals)
plt.show()

