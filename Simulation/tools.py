'''Thanks Niko, I bow to you.'''
import numpy as np
import matplotlib.pyplot as plt


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

def compPlot(tvec, solvec):
    ''' tvec is odesol.t
        solvec is the odesol.y'''

    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    fig1, axs = plt.subplots(3,1)

    axs[0].plot(tvec, x_vals)
    # axs[0].axis('equal')
    axs[0].set_title('X-Pos v. Time (s)', fontsize=10)
    axs[0].set_xlim(tvec[0], tvec[-1])
    # axs[0].set_ylim(np.min(x_vals), np.max(x_vals))

    axs[1].plot(tvec, y_vals)
    # axs[1].axis('equal')
    axs[1].set_title('Y-Pos v. Time (s)', fontsize=10)
    axs[1].set_xlim(tvec[0], tvec[-1])
    axs[1].set_ylim(np.min(y_vals), np.max(y_vals))

    axs[2].plot(tvec, z_vals)
    # axs[2].axis('equal')
    axs[2].set_title('Z-Pos v. Time (s)', fontsize=10)
    axs[2].set_xlim(tvec[0], tvec[-1])
    # axs[2].set_ylim(np.min(z_vals), np.max(z_vals))

    # return fig1

def Orbit3D(solvec, mu):
    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    n = np.linspace(0,2*np.pi,100)

    EarthX = -mu*np.cos(n)
    EarthY = -mu*np.sin(n)

    MoonX = (1-mu)*np.cos(n)
    MoonY = (1-mu)*np.sin(n)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(x_vals,y_vals,z_vals, c='r', s=.5, label='Craft')
    ax.scatter(mu, 0, 0, c='g', marker='x', s=20, label='Earth')
    ax.scatter(-(1-mu), 0, 0, c='b', marker='^', s=5, label='Moon')
    ax.scatter(0,0,0, c='m', marker='*')
    ax.plot3D(EarthX, EarthY, 0, c='g')
    ax.plot3D(MoonX, MoonY, 0, c='b')

    plt.axis('equal')
    ax.legend()
    plt.xlabel('X')
    plt.ylabel('Y')


def Orbit2D(solvec, mu):
    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    n = np.linspace(0,2*np.pi,100)

    EarthX = -mu*np.cos(n)
    EarthY = -mu*np.sin(n)

    MoonX = (1-mu)*np.cos(n)
    MoonY = (1-mu)*np.sin(n)

    fig = plt.figure()
    ax = plt.axes()
    ax.plot(x_vals,y_vals, c='r', label='Craft')
    ax.plot(-mu, 0, c='g', marker='x', label='Earth')
    ax.plot(1-mu, 0, c='b', marker='^', label='Moon')
    ax.plot(0,0, c='m', marker='*')
    ax.plot(EarthX, EarthY, c='g')
    ax.plot(MoonX, MoonY, c='b')

    plt.axis('equal')
    ax.legend()
    ax.set_xlim(-1,1)
    plt.xlabel('X')
    plt.ylabel('Y')