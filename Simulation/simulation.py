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
tbound = 18          # Final time
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
 
state0syn = [.715238, 0.19, .02, 0.06, .2, -.1 ] # IC in the synodic frame
# state0syn = [354494060.56376/c.lstar, 0.00001/c.lstar,  3880621.73670/c.lstar, 444.13436/ (c.lstar/c.tstar), -856.95527/ (c.lstar/c.tstar), -80.25872/ (c.lstar/c.tstar)]    #[ 0.7843, .3, 0, .2, -0.2, 0 ]

# state0syn = slib.StateDict('L3 Spider')

## Create the barycentric initial state vector
state0bar = ot.Syn2Bar(state0syn, t0, args={'Frame':'Barycentric'})

## Integrate initial state with a certian integrator and dynamics model/reference frame
SynSol = int.solve_ivp(fun=cr3bp.SynodicEOMs    , t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6)
# BarSol = int.solve_ivp(fun=cr3bp.BarycentricEOMs, t_span=[t0,tbound], y0=state0bar, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6)

# Move into Different Frame
# ECISol = ot.TransformSolFrame(BarSol, 'Bar2ECI')
BarSol = ot.SolTransform(SynSol, 'Syn2Bar')

## Plotting
pt.Orbit3D(SynSol.y, SynSol.t, c.mustar, args={'Frame':'Synodic'})     # plot 3d orbit in synodic frame
pt.Orbit3D(BarSol['y'], BarSol['t'], c.mustar, args={'Frame':'Barycentric'})     # plot 3d orbit in barycentric frame
# pt.Orbit3D(ECISol['y'], ECISol't'], c.mustar, args={'Frame':'Barycentric'}) # This actually is now in ECI

# pt.compPlot(BarSol['t'], BarSol['y'])

plt.show()