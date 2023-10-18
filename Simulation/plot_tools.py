"""Plotting tools file."""
import numpy as np
import matplotlib.pyplot as plt
import constants as c

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

def Orbit3D(solvec, time, mu, args={}):
    _args = {'Frame': 'Synodic'}
    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    traj = ax.scatter(x_vals,y_vals,z_vals, c=time, cmap = 'plasma', s=.5, label='Spacecraft')
    ax.scatter(0,0,0, c='m', marker='*', label='Barycenter')

    if args['Frame'] == 'Barycentric':
        n = np.linspace(0,2*np.pi,len(time))
        v = np.linspace(0, np.pi, len(time))

        re = c.earthD / c.lstar
        rm = c.moonD / c.lstar

        EarthX = -mu*np.cos(n)
        EarthY = -mu*np.sin(n)
        MoonX = (1-mu)*np.cos(n)
        MoonY = (1-mu)*np.sin(n)

        xe = re * np.outer(np.cos(n), np.sin(v)) - mu
        ye = re * np.outer(np.sin(n), np.sin(v)) + 0
        ze = re * np.outer(np.ones(np.size(n)), np.cos(v)) + 0

        xm = rm * np.outer(np.cos(n), np.sin(v)) + (1-mu)
        ym = rm * np.outer(np.sin(n), np.sin(v))
        zm = rm * np.outer(np.ones(np.size(n)), np.cos(v))
        # c=time, cmap = 'plasma',
        ax.scatter(EarthX, EarthY, 0, c='g', s=1)
        ax.scatter(MoonX, MoonY, 0, c='m', s=1)
        ax.plot_surface(xe,ye,ze)
        ax.plot_surface(xm,ym,zm)
        plt.title('Orbit in the Earth Centered Inertial Frame (ECRF).')
    
    else: 
        n = np.linspace(0,2*np.pi,100)
        v = np.linspace(0, np.pi, 100)

        re = c.earthD / c.lstar;
        rm = c.moonD / c.lstar;

        xe = re * np.outer(np.cos(n), np.sin(v)) + mu
        ye = re * np.outer(np.sin(n), np.sin(v)) + 0
        ze = re * np.outer(np.ones(np.size(n)), np.cos(v)) + 0

        xm = rm * np.outer(np.cos(n), np.sin(v)) - (1-mu)
        ym = rm * np.outer(np.sin(n), np.sin(v))
        zm = rm * np.outer(np.ones(np.size(n)), np.cos(v))

        ax.plot_surface(xe,ye,ze)
        ax.plot_surface(xm,ym,zm)
        plt.title('Orbit in the Earth-Moon Rotating Synodic Frame (EMR).')

    plt.axis('equal')
    ax.legend()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar(traj)

def Orbit2D(solvec, time, mu, args={}):
    _args = {'Frame':'Synodic'}
    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    fig = plt.figure()
    ax = plt.axes()
    traj = ax.scatter(x_vals,y_vals, c=time, cmap = 'plasma', s=.5, label='Spacecraft')
    ax.plot(-mu, 0, c='g', marker='x', label='Earth')
    ax.plot(1-mu, 0, c='b', marker='^', label='Moon')
    ax.plot(0,0, c='m', marker='*')


    # if args['Frame'] == 'Barycentric':
    #     n = np.linspace(0,2*np.pi,100)
    #     v = np.linspace(0, np.pi, 100)

    #     re = c.earthD / c.lstar
    #     rm = c.moonD / c.lstar

    #     EarthX = -mu*np.cos(n)
    #     EarthY = -mu*np.sin(n)
    #     MoonX = (1-mu)*np.cos(n)
    #     MoonY = (1-mu)*np.sin(n)

    #     xe = re * np.outer(np.cos(n), np.sin(v)) + mu
    #     ye = re * np.outer(np.sin(n), np.sin(v)) + 0

    #     xm = rm * np.outer(np.cos(n), np.sin(v)) - (1-mu)
    #     ym = rm * np.outer(np.sin(n), np.sin(v))

    #     ax.scatter(EarthX, EarthY, c='g')
    #     ax.scatter(MoonX, MoonY, c='b')
    #     ax.scatter(xe,ye)
    #     ax.scatter(xm,ym, c=time, cmap = 'plasma', label='Moon Orbit')

    plt.axis('equal')
    ax.legend()
    ax.set_xlim(-1,1)
    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.colorbar(traj)