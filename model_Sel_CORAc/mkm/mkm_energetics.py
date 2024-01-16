""" 
    functions for energy manipulation

"""

import numpy as np
from copy import deepcopy


def adjust_SHE_engs(engs, U_SHE, echem, alpha=0.5):
    """
      function to adjust (free) reaction energies accoring to CHE
 
      Parameters
      ----------
      engs : array/list of list
        energies of intermediates (can also be relative energies)
      U_SHE : float
        value of applied potential vs. SHE
      echem : list of booleans
        Boolean list to indicate which steps are electrochemical
      alpha : float
        symmetry factor for elementary reactions
 
      Returns
      -------
      engs : array/list of list
        energies of intermediates (can also be relative energies)
    
    """
    # apply CHE (RHE, SHE(ph=0))
    engs = deepcopy(engs)
    for i in range(len(engs)):
        if echem[i]:
            engs[i][1] += alpha * U_SHE
            engs[i][2] += U_SHE
    return(engs)


def make_rates(engs, pH, echem, PD_alk=True):
    """
      function to evaluate rates from (free) energies according to
      harmonic transition state theory with standard pre-factor
 
      Parameters
      ----------
      engs : array/list of list
        energies of intermediates (can also be relative energies)
      pH : float
        value of applied pH
      echem : list of booleans
        Boolean list to indicate which steps are electrochemical
      PD_alk : boolean
        if proton donor H2O (according to alkaline conditions)
 
      Returns
      -------
      rts : array 
        forward and backward rates for different reaction steps
    
    """
    # standard prefactor
    h = 4.135667e-15; kb = 8.6173331e-5; T = 300.
    nu = (kb * T) / h

    # make forward- and backward barriers
    dengs = np.array([[_dEf(eng), _dEb(eng)] for eng in engs]).flatten()
    rts = dengs * -1 / (kb*T); rts = np.exp(rts) * nu

    # include pH depenedency
    aH2O = 1.0; aOH = 10**(pH-14); aH3O = 10**(-1*pH)
    afw = aH2O; abw = aOH
    if not PD_alk:
        afw = aH3O; abw = aH2O

    for i in range(len(echem)):
        if echem[i]:
            rts[i*2] *= aH2O
            rts[i*2+1] *= aOH
    return rts


def _dEf(eng):
    """
      helper-function to compute forward reaction energy

    """
    return max(0.0, max(eng[1]-eng[0], eng[2]-eng[0]))


def _dEb(eng):
    """
      helper-function to compute backward reaction energy

    """
    return max(0.0, max(eng[1]-eng[2], eng[0]-eng[2]))


if __name__ == "__main__":
    pass

