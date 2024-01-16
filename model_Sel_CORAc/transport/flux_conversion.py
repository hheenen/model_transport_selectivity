""" 
    helper functions to interconvert flux/TOF
    parts of the transport code are based on mA/cm2 for 
    ease of working with output could be changed

"""

import numpy as np
from model_COR_Ac.data.constants import *


def fB(dE, T=300):
    """
      helper-function to compute rate constant based on harmonic TST

    """
    return(np.exp(-1.0 * max(0.0,dE) / (c_kB*T)))


def convert_TOF2flux(flux, roughness, facet=100):
    """
      helper-function to convert TOF into flux

    """
    # get facet dependent area per site
    AperSite = d_AperSite[facet]
    # flux input in tof --> need to be transformed to mol/cm2
    fluxM = flux/(AperSite*c_A2cm*c_A2cm * c_Na) * roughness # mol/cm2
    return(fluxM)


def f2j(f, nel, facet=100): # TOF --> mA/cm2
    """
      helper-function to convert: 
      input: TOF
      output: mA/cm2 (n_ECSA)

    """
    # get facet dependent area per site
    AperSite = d_AperSite[facet]
    
    # TOF --> mol/cm2
    fluxM = f/(AperSite*c_A2cm*c_A2cm * c_Na)# mol/cm2
    
    # mol/cm2 --> mA/cm2
    j = fluxM * c_Na * nel # to electrons
    j *= c_qe # to ampere
    j *= 1e3 # to mA
    return(j)
    

