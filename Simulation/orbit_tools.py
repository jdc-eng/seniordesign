"""Orbit calculations and reference frames tools."""
'''Thanks Niko, I bow to you.'''
import constants as c
import numpy as np

def invertKeplerTimeEquation(M,E0,e):
    """Invert the Kepler Time Equation
    INPUTS:
        M - mean anomaly
        E0 - initial guess for Eccentric Anomaly
        e - eccentricity
    OUTPUTS:
        E - Eccentric Anomaly
    """
    E=E0
    while (np.abs(M-E+e*np.sin(E))>1e-12):
        E = E - (M-E+e*np.sin(E))/(e*np.cos(E)-1)
    return E

def kepler2rv(a,e,Omega,I,omega,tmtp,mu):
    """ Convert Keplerian Elements to a position and velocity.
    INPUTS:
        a - semjax 
        e - eccentricity
        Omega - longitude of ascending node
        I - inclination
        omega - argument of periapsis
        tmtp - time since periapsis passage (t-tp)
        mu - gravitational parameter (units of this dictate units of output)
    OUTPUTS:
        r - inertial position vector
        v - inertial velocity vector
    """
    b = a*np.sqrt(1-e**2)
    Omega_rotation = np.array([[np.cos(Omega), np.sin(Omega), 0],[-np.sin(Omega), np.cos(Omega), 0],[0, 0, 1]])
    I_rotation = np.array([[1, 0, 0],[0, np.cos(I), np.sin(I)],[0, -np.sin(I), np.cos(I)]])
    omega_rotation = np.array([[np.cos(omega), np.sin(omega), 0],[-np.sin(omega), np.cos(omega), 0],[0, 0, 1]])
    pDCMi = np.matmul(omega_rotation,np.matmul(I_rotation,Omega_rotation))
    n = np.sqrt(mu/a**3)
    M = np.mod(n*(tmtp),2*np.pi)
    E = invertKeplerTimeEquation(M,M,e)
    r_p = np.transpose(np.array([a*(np.cos(E)-e),b*np.sin(E),0]))
    v_p = (a*n/np.linalg.norm(r_p))*np.array([-a*np.sin(E),b*np.cos(E),0])
    r = np.matmul(pDCMi.transpose(),r_p)
    v = np.matmul(pDCMi.transpose(),v_p)
    return r,v

def Syn2Bar(synstate, t, args={}):
    '''Assumes both states coincide at t=0'''
    x,y,z,vx,vy,vz = synstate
    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    # At = np.matrix([[np.cos(2*np.pi*t), -np.sin(2*np.pi*t), 0],
    #                 [np.sin(2*np.pi*t),  np.cos(2*np.pi*t), 0],
    #                 [0,                    0,               1]])

    _args = {'Frame': 'ECI'}
    if args['Frame'] == 'ECI':
        Vt = np.matrix([[vx - y],
                        [vy + x + c.mustar],
                        [  vz  ]])

        r = np.matrix([[x+c.mustar],
                    [y],
                    [z]])

        barstate = np.zeros(len(synstate))
        barstate[0:3] = np.matmul( At, synstate[0:3] )
        barstate[3], barstate[4], barstate[5]  = np.matmul( At, Vt)

    else:
        Vt = np.matrix([[vx - y],
                        [vy + x],
                        [  vz  ]])

        r = np.matrix([[x],
                    [y],
                    [z]])

        barstate = np.zeros(len(synstate))
        barstate[0:3] = np.matmul( At, r )
        barstate[3], barstate[4], barstate[5]  = np.matmul( At, Vt)
    return barstate


def ECI2Syn(state, t):
    rECI = np.matrix([[state[0]], [state[1]], [state[2]]])
    vECI = np.matrix([[state[3]], [state[4]], [state[5]]])
    
    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    J = np.matrix([[0,  1, 0],
                   [-1, 0, 0],
                   [0,  0, 0]])

    rEMR = np.transpose(np.matmul( np.transpose(At), rECI) + np.matrix([[c.mustar],[0],[0]]))
    vEMR = np.transpose(-np.matmul(np.transpose(At), vECI) - np.matmul(J, ( np.transpose(rEMR) + np.matrix([[c.mustar],[0],[0]]) )))
    return np.concatenate((rEMR, vEMR), axis = 1)
