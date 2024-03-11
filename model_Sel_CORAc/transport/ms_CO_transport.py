""" 
    Module with functions for "multiscale" simulations of 
    CO transport in combination with a microkinetic model
    --------------------------------------------------------------------------
    CO transport is modelled analytically following Fick's law of diffusion
    physically this is sufficient, no reaction of CO during transport is 
    expected

"""

import logging
from scipy.optimize import fmin
from model_Sel_CORAc.transport.flux_conversion import convert_TOF2flux
from model_Sel_CORAc.data.constants import *


def solve_MS_CO_analytical(mkin, U_SHE, pH, Lx, roughness):
    """
      function which solves a microkinetic model and CO diffusion
      according to Fick's 1st law in a multiscale Ansatz
      --> the transport limit's the flux of the reactant to the surface
      --> an analytical solution is equivalent to a numerical solution
          since no solution reactions are involved
 
      Parameters
      ----------
      mkin : mkm-object
        pre-created mkm object (see model_Sel_CORAc.mkm.model_surfrct)
      U_SHE : float
        value of applied potential vs. SHE
      pH : float
        value of applied pH
      Lx : float
        diffusion length in cm 
      roughness : float
        roughness of catalyst

      Returns
      -------
      ts : array 
        time evolution until steady-state
      ys : array
        coverages at timesteps (time, cov); y[-1,:] for steady-state
      pr : array
        production rate as defined through r_ind and p_ind
      p0 : float
        activity of CO which is proportional to its reduction in surface 
        concentration

    """
    y0 = [1.0] + [0.0]*mkin.n_int
    
    # minimize CO activity
    tol = 1e-6
    p0 = fmin(_cCO_diff, 1.0, args=(mkin, y0, U_SHE, pH, Lx, roughness), xtol=tol, ftol=tol, disp=False)#, maxfun=100)

    # check c0 convergence
    ts, ys, pr = mkin.run_kinetic_ODE(y0, U_SHE, pH, p0=[0.0,p0[0]]) #, 0.0])
    c0 = _dc_flux_CO(sum(pr), Lx, roughness)

    conv = abs(p0[0]/c0-1.0)
    if conv > 1e-4:
        logging.info("TRNSouterCO: Warning:  %.2f V, pH=%.2f, Lx=%.2e did not converge:\n"%(U_SHE,pH,Lx) + \
            "conv = %.3e; p0=%.3e vs. c0(analytic)=%.3e"%(conv, p0[0], c0))
   
    return(ts, ys, pr, p0[0])


def _cCO_diff(x, mkin, y0, U_SHE, pH, Lx, roughness):
    """
      helper-function for iterative optimization of flux and activity

    """
    # only surface kinetics relevant for total flux
    ts, ys, pr = mkin.run_kinetic_ODE(y0, U_SHE, pH, p0=[0.0,x[0]]) #, 0.0])
    flux = sum(pr)
    # match to CO diffusion following Fick's 1st law of diffusion
    c0 = _dc_flux_CO(flux, Lx, roughness)
    return(abs(c0 - x))


def _dc_flux_CO(flux, Lx, roughness):
    """
      helper-function matching a surface activity to
      a flux

    """
    # match to CO diffusion following Fick's 1st law of diffusion
    fluxM = convert_TOF2flux(flux, roughness) # in mol/cm2
    
    # units for CO conversion
    cb = 1.0/c_H_co # saturated CO=atm
    ccb = cb * 1e-3 # mol/cm3

    c0 = -1*((fluxM/(c_D_co*1e4))*Lx - ccb)

    # output in ml/cm3 --> std acitivity
    c0 *= 1e3
    a = c0 / cb
    return(a) # return activity



