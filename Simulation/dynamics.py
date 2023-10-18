"""CR3BP Dynamics File."""
import numpy as np
import constants as c

'''Sub-script of 1 refers to the larger primary, Earth. Sub-script 2 refers to the smaller primary, the Mooon.'''
def SynodicEOMs(t, state):
    x, y, z, vx, vy, vz = state
    mustar = c.mustar

    m1 = 1 - mustar
    m2 = mustar

    rb1 = -mustar
    rb2 = 1 - mustar

    r1 = np.sqrt( (x - rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2)**2 + y**2 + z**2 )        # Magnitude of r2-sat vector
    
    statedot = np.zeros(len(state))

    statedot[0:3] = state[3:]
    statedot[3] = x + 2*vy - m1* (x-rb1)/ r1**3 - m2* (x-rb2)/ r2**3
    statedot[4] = y - 2*vx - m1*    y   / r1**3 - m2*    y   / r2**3
    statedot[5] =           -m1*    z   / r1**3 - m2*    z   / r2**3

    return np.array(statedot)

def BarycentricEOMs(t, state):
    '''State vector is now 8x1; 7th element is theta=non-dim time'''
    x, y, z, vx, vy, vz = state
    mustar = c.mustar
    # mu1 = c.earthMass * c.G
    # mu2 = c.moonMass * c.G

    mu1 = 1 - mustar
    mu2 = mustar

    rb1 = - mustar
    rb2 = 1 - mustar

    r1 = np.sqrt( (x - rb1*np.cos(t))**2 + (y - rb1*np.sin(t))**2 + z**2 )      # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2*np.cos(t))**2 + (y - rb2*np.sin(t))**2 + z**2 )      # Magnitude of r2-sat vector
    
    statedot = np.zeros(len(state))

    statedot[0:3] = state[3:6]
    statedot[3] = -mu1*(x - rb1*np.cos(t))/r1**3 - mu2*(x + rb2*np.cos(t))/r2**3
    statedot[4] = -mu1*(y - rb1*np.sin(t))/r1**3 - mu2*(y - rb2*np.sin(t))/r2**3
    statedot[5] = -mu1*z/r1**3 - mu2*z/r2**3
    # statedot[6] = 1 # angular velocity of barycentric frame is non dimensional as 1 
    return np.array(statedot)


def two_body_ode( t, state, mu = c.EMmu ):
	# state = [ rx, ry, rz, vx, vy, vz ]

	r = state[ :3 ]
	a = -mu * r / np.linalg.norm( r ) ** 3

	return np.array( [
		state[ 3 ], state[ 4 ], state[ 5 ],
		    a[ 0 ],     a[ 1 ],     a[ 2 ] ] )