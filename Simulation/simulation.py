'''This file houses the simulation.'''
import constants as c
import orbit_tools as ot
import plot_tools as pt
import dynamics as cr3bp
from scipy import integrate as int
from scipy import fft
import numpy as np
import matplotlib.pyplot as plt
import state_library as slib

## Setup timestep of propagation
t0 = 0              # Initial tTime
tbound = 48         # Final time
tsteps = 8000
step = (tbound-t0)/tsteps

## Create initial state vector
# state0syn = [.15, -.25, -.06, .6, .5, .07 ] # IC in the synodic frame
# state0syn = [-0.85249293, 0, 0, 0, 0.33447503, -0.21472796]
state0syn = slib.StateDict('Soccer Ball')
# state0syn = ot.Bar2Syn(ot.ECI2Bar([-.1, .08, 0, -.2, -.1, -.08],0) ,0)

## Integrate initial state with a certian integrator and dynamics model/reference frame
SynSol = int.solve_ivp(fun=cr3bp.SynodicEOMs, t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6)
# sailSol = int.solve_ivp(fun=cr3bp.SailSynodicEOMs, t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6) 

## Plotting
pt.Orbit3D(SynSol.y, SynSol.t, c.mustar, args={'Frame':'Synodic'})                 # plot 3d orbit in synodic frame
# pt.compPlot(SynSol.t, SynSol.y)

ECRF = ot.SolTransform(ot.SolTransform({'y':SynSol.y, 't':SynSol.t}, 'Syn2Bar'), 'Bar2ECI')
pt.Orbit3D(ECRF['y'], SynSol.t, c.mustar, args={'Frame':'ECRF'})

# pt.Orbit3D(sailSol.y, sailSol.t, c.mustar, args={'Frame':'Synodic'})
# pt.compPlot(sailSol.t, sailSol.y)

# Tisser0 = ot.Tisserand(SynSol, c.G*c.earthMass)
# Tisser2 = ot.Tisserand(sailSol, c.G*c.earthMass)

# pt.genPlot(SynSol.t, Tisser0,  labels={'title':'Tisserand Criterion for Unperturbed Orbit','xlabel':'Time','ylabel':'Tisserand Criterion'}, legend={0:'Unperturbed'})
# pt.genPlot(sailSol.t, Tisser2, labels={'xlabel':'Time','ylabel':'Tisserand Criterion'}, legend={0:'Sail Acc. and Sun Gravity'})


sailAcc = np.zeros((3,len(SynSol.t)))
# for iter in range(len(SynSol.t)):
#     sailAcc[:,iter] = np.array([c.a_0*np.cos(c.ws*SynSol.t[iter]), -c.a_0*np.sin(c.ws*SynSol.t[iter]), 0])
Jacobi = ot.JacobiConstant(SynSol.y, sailAcc)
pt.JacobiPlot(Jacobi, SynSol.t, [-2.175,-2.05])
# print(ot.closeApproach(SynSol))

# pos = SynSol.y[0:3, :]
# pmags = np.zeros(np.size(pos))

# print(np.linalg.norm(pos[:,5002]))

plt.show()




# yf = fft.fft(Jacobi)
# xf = fft.fftfreq(tsteps, step)[:step//2]
# plt.plot(xf, 2.0/step * np.abs(yf[0:step//2]))

# plt.show()

