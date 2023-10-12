"""Store a few different types of ode solvers."""
import numpy as np
from scipy import integrate as int

def ode45(f,t,h):
    '''f is the function housing the dynamics equations you would like to integrate'''
    '''h = t[i+1]-t[i]'''
    f = np.zeros(len(t))
    f[0] = f0
    for i in range(0,len(t)-1):
        
        F1 = h*f(y[i],t[i])
        F2 = h*f((y[i]+F1/2),(t[i]+h/2))
        F3 = h*f((y[i]+F2/2),(t[i]+h/2))
        F4 = h*f((y[i]+F3),(t[i]+h))
    y[i+1] = y[i] + 1/6*(F1 + 2*F2 + 2*F3 + F4)
    return f

def rk8th(fun, t0, y0, t_bound, rtol=0.001, atol=1e-06, vectorized=False, first_step=None, **extraneous):

    odeSol = int.DOP853(fun, t0, y0, t_bound, rtol=0.001, atol=1e-06, vectorized=False, first_step=None, **extraneous)
    return odeSol


