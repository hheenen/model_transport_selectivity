""" 
    functions to describe flux in flow-channel following Sherwood relations
    this transport is only relevant for OH- in context of the acetate
    selectivity/ further liquid product transport can be ignored (no effect)

"""

import numpy as np
from model_COR_Ac.transport.flux_conversion import convert_TOF2flux
from model_COR_Ac.data.constants import *


def estimate_convective_pH_analytical(flx, pH_bulk, roughness, SR='plate'):
    """
      function which solves for the surface pH [OH- transport]
      (electrode surface / wall of flow channel) using Sherwood relations
 
      Parameters
      ----------
      flx : 2-D array
        TOF as given by mkm/transport module
      pH_bulk : float
        catholyte pH
      roughness : float
        roughness of catalyst

      Returns
      -------
      pH_surf : float
        pH at electrode surface / wall of flow channel

    """
    # get estimate for flux
    f_OH = convert_TOF2flux(flx[0]*4 + sum(flx[1:])*8, roughness)

    ve = 1 # ml/min see: 10.1038/s41929-018-0133-2 & 10.1038/s41929-019-0269-8
    ve *= 1/60. # cm3/s
    
    # now use Sherwood formula to estimate things
    # and get convective mass transfer coefficient ke | compute in cm
    # GDE dimensions see: 10.1038/s41929-018-0133-2 & 10.1038/s41929-019-0269-8
    A_cross = 0.5*0.15 # in cm
    L_cross = 0.15 # in cm
    Lc = 2.0 # in cm
    
    # viscosity and density of water
    mu_e = 8.9e-4 * 1e-2 # kg*m*s-2/m2 --> kg*cm*s-2/cm2
    rho_e = 0.997*1e-3 # kg/cm3
    # Schmidt
    Sc = mu_e/(rho_e * c_D_oh*1e4) # NOTE: unit cm2/s
    # Reynolds
    Re = L_cross * ve * rho_e/mu_e
    
    # (a) “absorption by a falling film of the dissolution of a solid wall into a falling film”
    #       from 10.26434/chemrxiv.13073348.v1, "Transport Phenomena", respectively
    if SR == 'film':
        Re_2 = Lc * ve * rho_e/mu_e
        ke = 1.017*(Lc/L_cross * Re_2 * Sc)**(1./3.) * (c_D_oh*1e4)/Lc
        dcOH = f_OH/ke 
        dcOH *= 1e3 # mol/l
        
        cOH = dcOH + 10**(pH_bulk-14.0) # total conc
        pH_surf = np.log10(cOH)+14.0
        ## print("a(2)", ke, pH_surf)

    # (b) Graetz-Nusselt-Problem --> dissolution in a tube --> VERY similar to above: (a) ke=0.001629646614630505; (b) 0.0013792762352202522
    #L_cross = L_cross * 2.0 # 
    elif SR == 'GNP':
        ke = 1.62 * (Re*Sc*L_cross/Lc)**(1./3.) * (c_D_oh*1e4)/L_cross
        dcOH = f_OH/ke 
        dcOH *= 1e3 # mol/l
        cOH = dcOH + 10**(pH_bulk-14.0) # total conc
        pH_surf = np.log10(cOH)+14.0
        ## print("b", ke, pH_surf)

    # (c) Plate -- Cussler book
    elif SR == 'plate':
        Re_2 = Lc * ve * rho_e/mu_e
        ke = 0.646 * (Re_2 * Sc)**(1./3.) * (c_D_oh*1e4)/Lc
        dcOH = f_OH/ke 
        dcOH *= 1e3 # mol/l
        cOH = dcOH + 10**(pH_bulk-14.0) # total conc
        pH_surf = np.log10(cOH)+14.0
        ## print("c", ke, pH_surf)
    
    return(pH_surf)


