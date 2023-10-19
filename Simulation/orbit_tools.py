"""Orbit calculations and reference frames tools. Mostly functions for frame transformations."""
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

def SolTransform(sol, rotation):
    _args = {'Bar2ECI': Bar2ECI,
            'Syn2Bar':  Syn2Bar,
            'Bar2Syn':  Bar2Syn
            }
    rotState = np.zeros( (6, len(sol.t)) )

    fun = _args[rotation]
    i = 0
    for t in sol.t:
        rotState[0,i], rotState[1,i], rotState[2,i], rotState[3,i], rotState[4,i], rotState[5,i]= fun(sol.y[:,i], t)
        i += 1
    RotSol = {'y':rotState, 't':sol.t}
    return RotSol

def Syn2Bar(synstate, t, args={}):
    '''Assumes both states coincide at t=0'''
    x,y,z,vx,vy,vz = synstate
    rSYN = np.matrix([[x],[y],[z]])
    vSYN = np.matrix([[vx],[vy],[vz]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    J = np.matrix([[0,  1, 0],
                   [-1, 0, 0],
                   [0,  0, 0]])

    rBAR = np.matmul( At, rSYN)
    vBAR = -np.matmul( At, np.matmul( J, rSYN ) ) + np.matmul( At, vSYN)

    barstate = np.zeros(len(synstate))
    barstate[0], barstate[1], barstate[2] = rBAR
    barstate[3], barstate[4], barstate[5] = vBAR
    return barstate

def Bar2Syn(barstate, t):
    x,y,z,vx,vy,vz = barstate
    rBAR = np.matrix([[x], [y], [z]])
    vBAR = np.matrix([[vx],[vy],[vz]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    J = np.matrix([[0,  1, 0],
                   [-1, 0, 0],
                   [0,  0, 0]])
    
    rEMR = np.matmul( np.transpose(At), rBAR )
    vEMR = np.matmul( np.transpose(At), ( vBAR + np.matmul( At, np.matmul( J, rEMR )) ) )

    synstate = np.zeros(len(barstate))
    synstate[0], synstate[1], synstate[2] = rEMR
    synstate[3], synstate[4], synstate[5] = vEMR
    return synstate

def Bar2ECI(state, t):
    mu = c.mustar
    rBAR = np.matrix([[state[0]], [state[1]], [state[2]]])
    vBAR = np.matrix([[state[3]], [state[4]], [state[5]]])

    rE   = np.matrix([[-c.mustar],       [0],        [0]])
    vE   = np.matrix([       [0],    [-mu*t],        [0]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    rECI = rBAR + np.matmul( At, rE )
    vECI = vBAR + np.matmul( At, vE )

    ECIstate = np.zeros(len(state))
    ECIstate[0], ECIstate[1], ECIstate[2] = rECI
    ECIstate[3], ECIstate[4], ECIstate[5] = vECI
    return ECIstate 

def ECI2Bar(state, t):
    mu = c.mustar
    rECI = np.matrix([[state[0]], [state[1]], [state[2]]])
    vECI = np.matrix([[state[3]], [state[4]], [state[5]]])

    rE   = np.matrix([[-c.mustar],       [0],        [0]])
    vE   = np.matrix([       [0],    [-mu*t],        [0]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    rBAR = rECI - np.matmul( At, rE )
    vBAR = vECI - np.matmul( At, vE )

    ECIstate = np.zeros(len(state))
    ECIstate[0], ECIstate[1], ECIstate[2] = rECI
    ECIstate[3], ECIstate[4], ECIstate[5] = vECI
    return ECIstate

# def Syn2ECI(synstate, t):
#     '''Assumes both states coincide at t=0'''
#     x,y,z,vx,vy,vz = synstate
#     At = np.matrix([[np.cos(t), -np.sin(t), 0],
#                     [np.sin(t),  np.cos(t), 0],
#                     [        0,          0, 1]])

#     Vt = np.matrix([[vx - y],
#                     [vy + x + c.mustar],
#                     [  vz  ]])

#     r = np.matrix([[x+c.mustar],
#                    [y],
#                    [z]])

#     barstate = np.zeros(len(synstate))
#     barstate[0:3] = np.matmul( At, synstate[0:3] )
#     barstate[3], barstate[4], barstate[5]  = np.matmul( At, Vt)


# def ECI2Syn(state, t):
#     rECI = np.matrix([[state[0]], [state[1]], [state[2]]])
#     vECI = np.matrix([[state[3]], [state[4]], [state[5]]])
#     rE   = np.matrix([[c.mustar],[0],[0]])

#     At = np.matrix([[np.cos(t), -np.sin(t), 0],
#                     [np.sin(t),  np.cos(t), 0],
#                     [        0,          0, 1]])
    
#     J = np.matrix([[0,  1, 0],
#                    [-1, 0, 0],
#                    [0,  0, 0]])
    
#     rEMR =  np.matmul( np.transpose(At), rECI ) + rE
#     vEMR = -np.matmul( np.transpose(At), vECI ) - np.matmul(J, ( rEMR + rE ))
#     synstate = np.zeros(len(state))
#     synstate[0], synstate[1], synstate[2] = rEMR
#     synstate[3], synstate[4], synstate[5] = vEMR
#     return synstate


