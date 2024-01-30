""" 
    Module with functions for "multiscale" simulations of 
    intermediate transport in combination with a microkinetic model

"""

import logging
import numpy as np
from scipy.optimize import fmin
from model_Sel_CORAc.transport.flux_conversion import convert_TOF2flux
from model_Sel_CORAc.data.constants import *
from model_Sel_CORAc.mkm.mkm_energetics import adjust_SHE_engs, make_rates
    


def solve_MS_X_analytical(mkin, U_SHE, pH, Dx, Lx, i_diss=1, i_rads=-1, roughness=1.0):
    """ 
      Function for "multiscale" simulations of analytical diffusion 
      in combination with a microkinetic model

      Parameters
      ----------
      mkin : mkm-object
        pre-created mkm object (see model_Sel_DRAR.mkm_variable.py)
      U_SHE : float
        value of applied potential vs. SHE
      pH : float
        value of applied pH
      Dx : float
        diffusion coefficient in m2/s
      Lx : float
        diffusion length in cm 
      i_diss = int
        indicator which pr rates go to transport related product
      i_rads = int
        indicate which rate is the re-adsorption
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

    pA = 1.0 # Dummy placeholder - activity of A never altered
    
    # initial range check
    c0max, increase_max = _search_min_diff(mkin, y0, U_SHE, pH, Dx, Lx, crange=1*10.**(np.arange(0,18)[::-1] * -1), 
        pA=pA, i_diss=i_diss, i_rads=i_rads, roughness=roughness) # check whole range always

    c0init = c0max[0] # Alternative: use c0init = 0 and tol = 1e-13 --> add warning if c0 ~ tol

    # use optimizer to solve for volatile product Cg near-surface concentration
    tol = c0init*1e-4 # tolerance two order of magnitude lower (percentage region)
    p0 = fmin(_c0_diff, c0init, args=(mkin, y0, U_SHE, pH, Dx, Lx, pA, i_diss, i_rads, roughness), xtol=tol, ftol=tol, disp=False)
    p0 = p0.tolist()+[pA] # pA pressure

    # check c0 convergence
    ts, ys, pr = _eval_flux_mkm(mkin, y0, U_SHE, pH, p0, i_diss, i_rads)
    c0 = compute_fick_c0(pr[0], Dx=Dx, Lx=Lx, pH=pH, roughness=roughness)
    
    conv = abs(p0[0]/c0-1.0)
    _warn_conv(conv, U_SHE, pH, Lx, p0, c0)

    # legacy due to bug in past: normalization of flux in case of leak 
    # in mass-conservation -- function left in to double check (prints warning)
    pr, rtot = _normalization_flux(mkin, y0, U_SHE, pH, pA, i_diss, i_rads, pr, conv)
    return(ts, ys, pr, p0[0]/rtot, p0)


def _c0_diff(x, mkin, y0, U, pH, Dx, Lx, pA, i_diss, i_rads, roughness):
    """
      helper-function solve mkm `mkin` for a surface concentration `x` and the
      analytical solution-reaction model for the resulting flux giving a surface
      concentration `c0`. The iterative adjustment of `x` and `c0` occurs via
      the fmin optimizer

    """
    # target function for fmin
    ts, ys, flux0 = _eval_flux_mkm(mkin, y0, U, pH, [x[0], pA], i_diss, i_rads)
    flux = max(0.0, flux0[0]) # for analytical solution non-negative flux is necessary
    c0 = compute_fick_c0(flux, Dx=Dx, Lx=Lx, pH=pH, roughness=roughness)

    return(abs(c0 - x))


#########################################################################
###   Routines for iterative solution (for analytical and numerical)  ###
#########################################################################

def _search_min_diff(mkin, y0, U, pH, Dx, Lx, crange, pA, i_diss, i_rads, conv_early=True, roughness=1.0):
    """
      helper-function screening a given range of concentration for a best fit;

    """
    cdiff = []; cref = []; cratio = []; pbreak = 0
    for x in crange:
        ts, ys, flux0 = _eval_flux_mkm(mkin, y0, U, pH, [x, pA], i_diss, i_rads)
        flux = max(0.0, flux0[0]) # for analytical solution non-negative flux is necessary
        c0 = compute_fick_c0(flux, Dx=Dx, Lx=Lx, pH=pH, roughness=roughness)
        cdiff.append(abs(c0-x))
        cref.append(c0)
        cratio.append(c0/x)
    imin = np.absolute(np.array(cratio) -1.0).argmin() # this is tested (same as cdiff)
    increase_max = np.all(np.array(cratio) > 1.0)
    ndiff = crange[imin-1:imin+2][[0,2]]
    return(ndiff, increase_max)


def _eval_flux_mkm(mkin, y0, U_SHE, pH, p0, i_diss=1, i_rads=-1):
    """
      helper-function to evaluate the flux from a microkinetic model object 
      "mkin" with adjustment of potential "U_SHE" and "pH" and activities "p0"; 
      the flux is re-evaluate since "mkm/ODEmkm" does not
      take gas-phase activity for production rates into account

    """
    # solve mkm
    ts, ys, pr = mkin.run_kinetic_ODE(y0, U_SHE, pH, p0=p0)

    # compute readsorption rate (U_SHE and pH independent)
    engsU = adjust_SHE_engs(mkin.engs0, U_SHE, mkin.echem)
    drt = make_rates(engsU, pH, [False]) 
    
    # evaluate flow (for one product hard-coded):
    #   i_diss > indicator which pr rates go to Cg-product
    #   p0 * rate_readsorption * empty sites
    flux0 = sum(pr[:i_diss]) - p0[0]*drt[i_rads]*ys[-1,0]
    # pr = [C_t, D_t]
    return(ts, ys, [flux0, *pr[i_diss:]])


def _normalization_flux(mkin, y0, U, pH, pA, i_diss, i_rads, pr, conv):
    """
      helper-function to normalize production rate
      this function is left in to double check, a missing mass-conservation
      was in this model due to an error in the mkm-equations;
      
      could in principle to depreciated (but doesn't hurt neither)

    """
    # post-normalize flux created from extra source term + Warning
    # (a) flux with NO extra outer pressure
    ts00, ys00, pr00 = _eval_flux_mkm(mkin, y0, U, pH, [0.0, pA], i_diss, i_rads)
    # (b) add Warning/message if outer flux > 0.01%
    rtot = sum(pr) / sum(pr00)
    if rtot > 1.0001 and rtot < 10.: # left in to double check
        print("TRNSouter: Warning: increase in total rate through" + \
            "transport %.2e vs. %.2e (ratio %.3f)"%(sum(pr), sum(pr00), rtot-1.0))
    elif rtot > 10: # total rate increase > 10 counted as not converged
        print("TRNSouter: Failed convergence: increase in total rate through" + \
            "transport %.2e vs. %.2e (ratio %.3f)"%(sum(pr), sum(pr00), rtot-1.0))
    # (c) renormalize pr (not ys -- although should)
    pr /= rtot
    # (d) give np.nan to non-converged values
    if conv > 0.5 or rtot > 10:
        pr[:] = np.nan
    return(pr, rtot)


def _warn_conv(conv, U, pH, Lx, p0, c0):
    if conv > 1e-4:
        print("TRNSouter: Warning:  %.2f V, pH=%.2f, Lx=%.2e did not converge:\n"%(U,pH,Lx) + \
            "conv = %.3e; p0=%.3e vs. c0(analytic)=%.3e"%(conv, p0[0], c0))


def compute_fick_c0(flux, Dx, Lx, pH, roughness=1.0): # steady-state solution Fick
    """
      function for analytical nearsurface concentration

    """
    # conversions in cm2/s
    fluxM = convert_TOF2flux(flux, roughness)

    # note on conversions 1e4 (m2/s --> cm2/s); 1e3 (mol/m3 * s --> mol/l * s)
    c0a = get_analytical_conc0(fluxM, Lx, Dx*1e4)
    
    # output in ml/cm3 = 1e3 * std acitivity...
    c0a *= 1e3
    return(c0a)


def get_analytical_conc0(cgrd, Lx, D, cbulk=0):
    """
      helper-function for analytical near-surface concentration of Cg in
      diffusion-only case

    """
    ## just Fick's 1st law of diffuison: cgrd = D*dy/dx
    dy = (cgrd*Lx)/D
    y0 = dy + cbulk
    return(y0)




