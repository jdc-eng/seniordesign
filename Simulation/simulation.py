'''This file houses the simulation.'''
import constants as c
import orbit_tools as ot
import plot_tools as pt
import dynamics as cr3bp
from scipy import integrate as int
import numpy as np
import matplotlib.pyplot as plt
import state_library as slib

## Setup timestep of propagation
t0 = 0              # Initial tTime
tbound = 12          # Final time
tsteps = 4000
step = (tbound-t0)/tsteps

# Create state vector
# r, v = ot.kepler2rv(
#         a     = 100000 * 1000,
#         e     = 0.6,
#         Omega =  79,
#         I     =   0,
#         omega =  70,
#         tmtp  =   0,
#         mu    = c.G*c.earthMass)

# state0ECI = np.concatenate( ( r / c.lstar,  v / (c.lstar/c.tstar) ) )
# state0syn = np.asarray(ot.ECI2Syn(state0ECI, t0))
 
# state0syn = [ -.49, .7, -.3, .2,-0.4, 0.4 ] # IC in the synodic frame
# state0syn = [ 0.7843, .3, 0, .2, -0.2, 0 ]
# state0bar = ot.Syn2Bar(state0syn, t0, args={'Frame':'ECI'})

state0syn = slib.StateDict('Turbine')
state0bar = ot.Syn2Bar(state0syn, t0, args={'Frame':'Barycentric'})


## Integrate initial state with a certian integrator and dynamics model

SynSol = int.solve_ivp(fun=cr3bp.SynodicEOMs    , t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6)
BarSol = int.solve_ivp(fun=cr3bp.BarycentricEOMs, t_span=[t0,tbound], y0=state0bar, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6)

## Move into ECI Frame
# BarSol = np.zeros( (7, len(SynSol.t)) )
# i = 0
# for t in SynSol.t:
#     BarSol[0,i], BarSol[1,i], BarSol[2,i], BarSol[3,i], BarSol[4,i], BarSol[5,i]= ot.Syn2Bar(SynSol.y[:,i], t, args={'Frame':'ECI'})
#     BarSol[6,i] = t
#     i += 1


## Plotting

pt.Orbit3D(SynSol.y, SynSol.t, c.mustar, args={'Frame':'Synodic'})     # plot 3d orbit in synodic frame
pt.Orbit3D(BarSol.y, BarSol.t, c.mustar, args={'Frame':'Barycentric'})     # plot 3d orbit in barycentric frame
# pt.Orbit3D(BarSol[0:5,:], BarSol[6,:], c.mustar, args={'Frame':'Barycentric'}) # This actually is now in ECI
# pt.Orbit2D(BarSol[0:5,:], BarSol[6,:], c.mustar, args={'Frame':'Barycentric'})

plt.show()