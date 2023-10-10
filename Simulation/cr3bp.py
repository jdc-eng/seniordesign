"""CR3BP Dynamics File."""
import numpy as np
import constants as c

# put sail's force/acceleration into its 3 components, add them into the eom for 3bp

'''Sub-script of 1 refers to the larger primary, Earth. Sub-script 2 refers to the smaller primary, the Mooon.'''

# pass in a state vector of position and velocity vector components

def SynodicEOMs(t, state):
    x, y, z, vx, vy, vz = state
    mustar = c.mustar
    
    r1 = np.sqrt((x-mustar)**2+y**2+z**2)          # Magnitude of r1-sat vector
    r2 = np.sqrt((x-mustar+1)**2+y**2+z**2)        # Magnitude of r2-sat vector

    statedot = np.zeros(len(state))

    statedot[0:3] = state[3:]
    statedot[3] = x + 2*vy - (1-mustar)*(x-mustar)/r1**3-mustar*(x+1-mustar)/r2**3
    statedot[4] = y-2*vx-(1-mustar)*y/r1**3 - mustar*y/r2**3
    statedot[5] = -(1-mustar)*z/r1**3-mustar*z/r2**3
    return np.array(statedot)

