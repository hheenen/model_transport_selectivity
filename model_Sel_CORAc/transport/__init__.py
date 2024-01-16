"""
This module contains all functions to handle the transport and integrate it
with the kinetic model from model_COR_Ac.mkm in a multi-scale ansatz

most functionality is summarized in the wrapper function `solve_MS_model`

"""

import numpy as np

from model_COR_Ac.transport.ms_CO_transport import solve_MS_CO_analytical
from model_COR_Ac.transport.ms_Ac_transport import solve_MS_ketene_analytical, \
                                                   solve_MS_ketene_numerical
from model_COR_Ac.transport.convection_analytical import estimate_convective_pH_analytical


def solve_MS_model(mkin, U_SHE, pH, Lx, diffmode, i_diss=1, i_rads=-1, roughness=1.0, tol_diff=None):
    """ 
      Wrapper-function for "multiscale" simulations of CO, ketene and OH- 
      transport limitations in combination with a microkinetic model.
      analytical and numerical solvers can be used

      Parameters
      ----------
      mkin : mkm-object
        pre-created mkm object (see model_COR_Ac.mkm.model_surfrct)
      U_SHE : float
        value of applied potential vs. SHE
      pH : float
        value of applied pH
      Lx : float
        diffusion length in cm 
      diffmode : str
        keyword for mode of ketene transport (diff / diff-rct / diff-rct-num)
      pCO : float
        activity of CO at catalyst surface
      i_diss = int
        indicator which pr rates go to transport related product
      i_rads = int
        indicate which rate is the re-adsorption
      roughness : float
        roughness of catalyst
      tol_diff : float
        tolerance for numerical solver to converge

      Returns
      -------
      ts : array 
        time evolution until steady-state
      ys : array
        coverages at timesteps (time, cov); y[-1,:] for steady-state
      pr : array
        production rate as defined through r_ind and p_ind
      p0 : list
        activity / concentration of [ketene, CO]
      Cx : dictionary
        concentration profile of ketene and OH

      --------------------------------------------------------------------------
      The solution reaction and diffusion is modelled using a finite difference
      scheme (see transport/diffrct_numerical_1D2cmp.py) the implementation is 
      not very efficient (equidistant grid) and can only be used for small 
      diffusion lengths which, however, suffice in this work

    """
    
    # function to solve both CO and ketene surface concentration analytically or numerically

    # Lx can be given individually (gas and solution side)
    if type(Lx) != list:
        Lx = [Lx, Lx]
    
    # (1) solve CO activity - mass transport limitations of 
    ts, ys, pr, pCO = solve_MS_CO_analytical(mkin, U_SHE, pH, Lx[0], roughness)
    
    # (2) solve OH concentration due to convection at electrode
    pH_local = estimate_convective_pH_analytical(pr, pH, roughness)
    if diffmode != 'diff-rct-num':
        # (3) analytical solution for ketene solution reaction
        ts, ys, pr, pKt, praw = solve_MS_ketene_analytical(mkin, U_SHE, pH_local, 
            Lx[1], diffmode, pCO, i_diss, i_rads, roughness)
        Cx = None
    else:
        # (3) numerical solution for ketene solution reaction
        ts, ys, pr, pKt, praw, Cx = solve_MS_ketene_numerical(mkin, U_SHE, pH_local, 
            Lx[1], diffmode, pCO, i_diss, i_rads, roughness, 
            pklfile='dat_bk_num_flux.pkl', tol_diff=tol_diff)
    
    # pr here in TOF, [pKt, pCO] in total, need to scale TOF with roughness | for consistend output
    pr = (np.array(pr)*roughness).tolist()

    return(ts, ys, pr, [pKt, pCO], praw, Cx)

