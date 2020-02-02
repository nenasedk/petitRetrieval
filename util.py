import numpy as np
import math as math
from sys import platform
import os
import threading, subprocess
from scipy.special import erfcinv

SQRT2 = math.sqrt(2.)

# Sanity checks on parameter ranges
def b_range(x, b):
    if x > b:
        return -np.inf
    else:
        return 0.

def a_b_range(x, a, b):
    if x < a:
        return -np.inf
    elif x > b:
        return -np.inf
    else:
        return 0.

# Convert from emission flux to measured flux at earth
def Surf_To_Meas(flux,p_rad,dist):
    # Converts ergs cm^-2 s^-1 Hz^-1 to uJy 
    # Computes from surface emission flux to measured flux at distance dist
    # Distances must be in the same units
    m_flux = flux * p_rad**2/dist**2
    return m_flux*1e29

# Priors stolen from https://github.com/JohannesBuchner/MultiNest/blob/master/src/priors.f90
def log_prior(cube,lx1,lx2):
    return 10**(lx1+cube*(lx2-lx1))

def uniform_prior(cube,x1,x2):
    return x1+cube*(x2-x1)

def gaussian_prior(cube,mu,sigma):
    return mu + sigma*SQRT2*erfcinv(2.0*(1.0 - cube))
    #return -(((cube-mu)/sigma)**2.)/2.

def log_gaussian_prior(cube,mu,sigma):
    bracket = sigma*sigma + sigma*SQRT2*erfcinv(2.0*cube)
    return bracket

def show(filepath): 
    """ open the output (pdf) file for the user """
    if os.name == 'mac' or platform == 'darwin': subprocess.call(('open', filepath))
    elif os.name == 'nt' or platform == 'win32': os.startfile(filepath)
    elif platform.startswith('linux') : subprocess.call(('xdg-open', filepath))
