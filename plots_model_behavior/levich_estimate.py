#!/usr/bin/env python

"""

This script only produces 1D cuts for each parameter

"""

import os, sys
import numpy as np
from copy import deepcopy


if __name__ == "__main__":
    
    #####################################################
    #### one dimensional diffusion problem -- Levich ####
    #####################################################
    
    # estimate d_levic --> for comparison of tested diffusion layer thicknesses
    D_CO = 2.03e-5 * 1e-4 # cm2/s --> m2/s
    D_range = np.array([1.0, 2.2]) * 1e-5 * 1e-4 # 
    # dynamic --> kinematic viscosities (mPa s) / (g/cm3) = m2/s *1e-6
    visc =  (0.9980/1.04652 +  1.1200/1.091221)/2.0 * 1e-6 # mPa s --> Ns/m2 for KOH 10.1021/je000019h
    v_range = np.array([0.5, 1.5]) * 1e-6 # --> temperature range of water

    rpms = [100, 400, 800, 1600, 2500]
    ds = []
    for rpm in rpms:
        for D in D_range:
            for v in v_range:
                o = rpm * np.pi*2/60
                dlevic = 1.6126 * D_CO**(1./3.) * visc**(1./6.) * o**(-1./2.) *1e2 # in cm
                ds.append(dlevic * 1e-2 *1e6)

    print(min(ds), max(ds), 'in mu')
