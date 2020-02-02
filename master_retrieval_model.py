##----------------------------------------
## Master retrieval model
##----------------------------------------

import numpy as np
import sys
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc

def calc_MMW(abundances):

    MMWs = {}
    MMWs['H2'] = 2.
    MMWs['He'] = 4.
    MMWs['H2O'] = 18.
    MMWs['CH4'] = 16.
    MMWs['CO2'] = 44.
    MMWs['CO'] = 28.
    MMWs['Na'] = 23.
    MMWs['K'] = 39.
    MMWs['NH3'] = 17.
    MMWs['HCN'] = 27.
    MMWs['C2H2'] = 26.
    MMWs['PH3'] = 34.
    MMWs['H2S'] = 34.
    MMWs['VO'] = 67.
    MMWs['TiO'] = 64.
    MMWs['H2S'] = 34.

    MMW = 0.
    for key in abundances.keys():
        if key == 'CO_all_iso':
            MMW += abundances[key]/MMWs['CO']
        else:
            MMW += abundances[key]/MMWs[key]
    
    return 1./MMW
    
####################################################################################
####################################################################################
####################################################################################

def retrieval_model_plain(rt_object, temperature_parameters, log_g, log_P0, \
                              R_pl, ab_metals):

    gravity = 1e1**log_g    
    
    # Create temperature model
    press, temp = nc.make_press_temp(temperature_parameters) # pressures from low to high

    abundances = {}
    metal_sum = 0.
    for name in ab_metals.keys():
        abundances[name] = np.ones_like(press)*1e1**ab_metals[name]
        metal_sum += 1e1**ab_metals[name]

    abH2He = 1. - metal_sum
    abundances['H2'] = np.ones_like(press)*abH2He*0.75
    abundances['He'] = np.ones_like(press)*abH2He*0.25
            
    MMW = calc_MMW(abundances)
        
    rt_object.calc_flux(temp, abundances, gravity, MMW)
        
    return nc.c/rt_object.freq, rt_object.flux
