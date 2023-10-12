'''This file houses the simulation.'''
import constants as c
import tools
import cr3bp
from scipy import integrate as int
import numpy as np
import matplotlib.pyplot as plt
import state_library as slib

## Setup timestep of propagation
t0 = 0              # Initial tTime
tbound = 8    # Final time
tsteps = 10000
step = (tbound-t0)/tsteps

## Create state vector
r, v = tools.kepler2rv(
        a     = 250000 * 1000,
        e     = .05,
        Omega = 200,
        I     =  40,
        omega = 30,
        tmtp  =   0,
        mu    = c.G*c.earthMass)

state = np.concatenate( ( r / c.lstar,  v / (c.lstar/c.tstar) ) )
# state = [0.5, 0, 0, 0, -0.1, 0] # Pointy star
# state = slib.StateDict('Pointy Star')



## Integrate those dynamics with a ce[0rtian solver and dynamics model
odesol = int.solve_ivp(fun=cr3bp.SynodicEOMs, t_span=[t0,tbound], y0=state, method='DOP853', max_step = step, atol=1e-12, rtol=1e-9)


## Plotting
# tools.compPlot(odesol.t, odesol.y)      # plot each position component v time
tools.Orbit3D(odesol.y, c.mustar)       # plot 3d orbit in synodic frame
# tools.Orbit2D(odesol.y, c.mustar)

plt.show()