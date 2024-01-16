""" 
    Module with functions for "multiscale" simulations of numerical ketene 
    transport and solution reaction in combination with a microkinetic model
    --------------------------------------------------------------------------
    The solution reaction and diffusion is modelled using a finite difference
    scheme - the implementation is not very efficient (equidistant grid)
    and can only be used for small diffusion lengths (enough in this example)

"""

import numpy as np
import os, logging
from scipy.optimize import fmin
from model_COR_Ac.transport.flux_conversion import fB, convert_TOF2flux
from model_COR_Ac.transport.diffrct_numerical_1D2cmp import diffrct1d2cmp
from model_COR_Ac.tools import load_pickle_file, write_pickle_file, _make_hash
from model_COR_Ac.mkm.mkm_energetics import adjust_SHE_engs, make_rates
from model_COR_Ac.data.constants import *


#########################################################################
#####               function for analytical solution                #####
#########################################################################

def solve_MS_ketene_numerical(mkin, U_SHE, pH, Lx, diffmode, pCO, i_diss, 
    i_rads, roughness, pklfile=None, tol_diff=None):
    """ 
      Function for "multiscale" simulations of numerical ketene transport and 
      solution reaction in combination with a microkinetic model

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
        keyword for mode of ketene transport (diff / diff-rct)
      pCO : float
        activity of CO at catalyst surface
      i_diss = int
        indicator which pr rates go to transport related product
      i_rads = int
        indicate which rate is the re-adsorption
      roughness : float
        roughness of catalyst
      pklfile : str
        name of pklfile used to store converged results of transport solver
        (meant to speed-up - but not confirmed)
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
      p0 : float
        activity of CO which is proportional to its reduction in surface 
        concentration

      --------------------------------------------------------------------------
      The solution reaction and diffusion is modelled using a finite difference
      scheme (see transport/diffrct_numerical_1D2cmp.py) the implementation is 
      not very efficient (equidistant grid) and can only be used for small 
      diffusion lengths which, however, suffice in this work

    """
    y0 = [1.0] + [0.0]*mkin.n_int
    
    # initial range check: determine search window and tolerance
    c0max, increase_max = _search_min_diff(mkin, y0, U_SHE, pH, Lx, 
        crange=1*10.**(np.arange(0,18)[::-1] * -1), pCO=pCO, diffmode='diff-rct', 
        i_diss=i_diss, i_rads=i_rads, roughness=roughness)
    c0init = min(1e-3, c0max[0]) # 1e-3 mol/cm3 (1mol/l) as maximum range
    tol = c0init*1e-4 # accuracy = 0.01 %
    logging.debug("TRNSouter: tolerance for convergence fmin %.4e"%tol)

    # initital guess for flux
    ts, ys, pr_guess = _eval_flux_mkm(mkin, y0, U_SHE, pH, [c0init, pCO], 
        i_diss, i_rads)
    # load previous solutions if exist
    Cinit = _recover_saved_solutions(pklfile, Lx, [c_D_kt, 2*c_D_oh], c_k_sol,
        _make_hash(np.array(mkin.engs0)), pH, pr_guess)

    # initialize finite difference solver with system specific properties
    drnum = diffrct1d2cmp([c_D_kt, 2*c_D_oh], [c_k_sol, c_k_sol], Lx=Lx*1e-2, 
        dt=None, dx=5e-8, log=True, tol=tol_diff)
    drnum.Cx = Cinit # set initial guess

    p0 = fmin(_c0_diff_num, c0init, args=(mkin, y0, U_SHE, pH, drnum, pCO, 
        i_diss, i_rads, roughness), xtol=tol, ftol=tol, disp=False)#, maxfun=100)
    
    p0 = p0.tolist()+[pCO] #CO pressure

    # save last solution for initial guesses
    _save_num_solution(pklfile, Lx, [c_D_kt, 2*c_D_oh], c_k_sol, 
        _make_hash(np.array(mkin.engs0)), pH, pr_guess, drnum.Cx)

    # check c0 convergence
    ts, ys, pr = _eval_flux_mkm(mkin, y0, U_SHE, pH, p0, i_diss, i_rads)
    drnum.Cx = _recover_saved_solutions(pklfile, Lx, [c_D_kt, 2*c_D_oh], c_k_sol, 
        _make_hash(np.array(mkin.engs0)), pH, pr)
    c0 = _eval_diffrct_num(drnum, pr, pH, roughness)

    conv = abs(p0[0]/c0-1.0)
    _warn_conv(conv, U_SHE, pH, Lx, diffmode, p0, c0)
    
    # normalization due to iterative procedure -- mass-conservation leak
    # in the past -- function left in to double check (print's warning)
    pr, rtot = _normalization_flux(mkin, y0, U_SHE, pH, pCO, i_diss, i_rads, pr, conv)
    
    # nullify non-converged solutions
    if np.all(np.isnan(pr)):
        p0 = [np.nan]*len(p0)
        drnum.Cx[:,:] = np.nan

    return(ts, ys, pr, p0[0]/rtot, p0, drnum.Cx)


def _c0_diff_num(x, mkin, y0, U, pH, drnum, pCO, i_diss, i_rads, roughness):
    """
      helper-function solve mkm `mkin` for a surface concentration `x` and the
      transport model `drnum` for the resulting flux `flux0` giving a surface
      concentration `c0`. The iterative adjustment of `x` and `c0` occurs via
      the fmin optimizer

    """
    # target function for fmin
    ts, ys, flux0 = _eval_flux_mkm(mkin, y0, U, pH, [x[0], pCO], i_diss, i_rads)

    c0 = _eval_diffrct_num(drnum, flux0, pH, roughness)
    if np.isnan(c0): #failure in integration (timestep)
        raise Exception("integration failed to to timestep, bad guess?")
        return(0.0)
    
    logging.debug("TRNSinner: convergence fmin %.4e"%abs(c0 - x))

    return(abs(c0 - x))


def _eval_diffrct_num(drnum, flux0, pH, roughness):
    """
      helper-function to solve numerical model `drnum` for specific flux (TOF)

    """
    # add total flux for OH
    # NOTE: conversions in cm2/s
    fluxAc = convert_TOF2flux(flux0[0], roughness)
    fluxOH = convert_TOF2flux(flux0[0]*4 + sum(flux0[1:])*8, roughness)
    
    # solve numerical model
    cOH = 10**(pH-14) * 1e3 # mol/m3
    drnum.converge_num([fluxAc*1e4, fluxOH*1e4], [0.0, cOH], initC=drnum.Cx, stout=False)
    dat = drnum.get_concentration_profile()
    c0 = dat['conc'][0,0] * 1e-3 # into mol/l --> activity

    return(c0)


def _recover_saved_solutions(pklfile, Lx, D, kexp, ehsh, pH, flux):
    """
      helper-function to retrieve numerical concentration profile from pkl-file

    """
    if pklfile == None or not os.path.isfile(pklfile):
        return(np.array([]))
    else:
        d = load_pickle_file(pklfile)
        idkey = "%.2e_%.2e_%.2e_%.2e_%s_%.3f"%(Lx, *D, kexp, ehsh, pH)
        if idkey in d:
            flist = [[float(ff)for ff in f.split('_')] for f in d[idkey].keys()]
            flist = np.array(flist) # now find closest match
            ind = np.absolute(flist - flux).sum(axis=1).argmin()
            fkey = flist[ind]
            return(d[idkey]['%.2e_%.2e'%tuple(fkey)])
        else:
            return(np.array([]))


def _save_num_solution(pklfile, Lx, D, kexp, ehsh, pH, flux, Cx):
    """
      helper-function to save numerical concentration profile to pkl-file

    """
    if pklfile != None:
        idkey = "%.2e_%.2e_%.2e_%.2e_%s_%.3f"%(Lx, *D, kexp, ehsh, pH)
        d = {}
        if os.path.isfile(pklfile):
            d = load_pickle_file(pklfile)
        if idkey not in d:
            d.update({idkey:{}})
        d[idkey].update({'%.2e_%.2e'%tuple(flux):Cx})
        write_pickle_file(pklfile,d)
        

#########################################################################
#####               function for analytical solution                #####
#########################################################################

def solve_MS_ketene_analytical(mkin, U_SHE, pH, Lx, diffmode, pCO=1.0, i_diss=1, i_rads=-1, roughness=1.0):
    """ 
      Function for "multiscale" simulations of analytical ketene transport and 
      solution reaction in combination with a microkinetic model

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
        keyword for mode of ketene transport (diff / diff-rct)
      pCO : float
        activity of CO at catalyst surface
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


      This analytical solution is only an alternative to the more accurate
      numerical solution (results in the manuscript are based on the numerical
      version)
      --------------------------------------------------------------------------
      The solution reaction and diffusion is modelled using an analytical 
      derivation of the reaction-diffusion system. Here, a constant OH- 
      concentration is assumed. This assumption is justified since 
      [OH-] >> [ketene] (see manuscript) and therefore results differ only
      marginally to numerical solution

    """

    y0 = [1.0] + [0.0]*mkin.n_int
    
    # initial range check
    c0max, increase_max = _search_min_diff(mkin, y0, U_SHE, pH, Lx, crange=1*10.**(np.arange(0,18)[::-1] * -1), 
        pCO=pCO, diffmode=diffmode, i_diss=i_diss, i_rads=i_rads, roughness=roughness) # check whole range always
    c0init = c0max[0]

    # us optimizer to solve for ketene surface concentration
    tol = c0init*1e-4 # tolerance two order of magnitude lower (percentage region)
    p0 = fmin(_c0_diff, c0init, args=(mkin, y0, U_SHE, pH, Lx, diffmode, pCO, i_diss, i_rads, roughness), xtol=tol, ftol=tol, disp=False)#, maxfun=100)
    p0 = p0.tolist()+[pCO] #CO pressure

    # check c0 convergence
    ts, ys, pr = _eval_flux_mkm(mkin, y0, U_SHE, pH, p0, i_diss, i_rads)
    c0s = compute_ketene_c0(pr[0], Lx=Lx, pH=pH, roughness=roughness)
    c0 = c0s[diffmode]
    
    conv = abs(p0[0]/c0-1.0)
    _warn_conv(conv, U_SHE, pH, Lx, diffmode, p0, c0)

    # normalization due to iterative procedure -- mass-conservation leak
    # in the past -- function left in to double check (print's warning)
    pr, rtot = _normalization_flux(mkin, y0, U_SHE, pH, pCO, i_diss, i_rads, pr, conv)
    return(ts, ys, pr, p0[0]/rtot, p0)


def _c0_diff(x, mkin, y0, U, pH, Lx, diffmode, pCO, i_diss, i_rads, roughness):
    """
      helper-function solve mkm `mkin` for a surface concentration `x` and the
      analytical solution-reaction model for the resulting flux giving a surface
      concentration `c0`. The iterative adjustment of `x` and `c0` occurs via
      the fmin optimizer

    """
    # target function for fmin
    ts, ys, flux0 = _eval_flux_mkm(mkin, y0, U, pH, [x[0], pCO], i_diss, i_rads)
    flux = max(0.0, flux0[0]) # for analytical solution non-negative flux is necessary
    c0s = compute_ketene_c0(flux, Lx=Lx, pH=pH, roughness=roughness)
    c0 = c0s[diffmode]
    return(abs(c0 - x))


#########################################################################
###   Routines for iterative solution (for analytical and numerical)  ###
#########################################################################

def _search_min_diff(mkin, y0, U, pH, Lx, crange, pCO, diffmode, i_diss, i_rads, conv_early=True, roughness=1.0):
    """
      helper-function screening a given range of concentration for a best fit;
      based on analytical solution-reaction ("compute_ketene_c0") as an 
      approximation; 

    """
    cdiff = []; cref = []; cratio = []; pbreak = 0
    for x in crange:
        ts, ys, flux0 = _eval_flux_mkm(mkin, y0, U, pH, [x, pCO], i_diss, i_rads)
        flux = max(0.0, flux0[0]) # for analytical solution non-negative flux is necessary
        c0s = compute_ketene_c0(flux, Lx=Lx, pH=pH, roughness=roughness)
        c0 = c0s[diffmode]
        cdiff.append(abs(c0-x))
        cref.append(c0)
        cratio.append(c0/x)
    imin = np.absolute(np.array(cratio) -1.0).argmin() # this is tested (same as cdiff)
    increase_max = np.all(np.array(cratio) > 1.0)
    if imin == 0: # TODO if these don't occur anymore (since usage of fmin --> delete) ! wait for boundary condition test rationalization
        logging.debug("TRNSrange: input range out of range low-end")
        ndiff = [crange[0]-0.5*crange[0],crange[1]]
    elif imin == crange.size-1 or increase_max:
        logging.debug("TRNSrange: input range out of range high-end")
        ndiff = [crange[-1], crange[-1]+0.5*crange[0]]
    # encompassing crange (around minimum)
    else: # standard solution
        ndiff = crange[imin-1:imin+2][[0,2]]
    return(ndiff, increase_max)


def _eval_flux_mkm(mkin, y0, U_SHE, pH, p0, i_diss=1, i_rads=-1):
    """
      helper-function to evaluate the flux from a microkinetic model object 
      "mkin" with adjustment of potential "U_SHE" and "pH" and CO and ketene
      activities "p0"; the flux is re-evaluate since in "mkm/ODEmkm" does not
      take gas-phase activity for production rates into account

    """
    # solve mkm
    ts, ys, pr = mkin.run_kinetic_ODE(y0, U_SHE, pH, p0=p0)
    # compute flux - normally should be in production rate (no gas-conc internally)

    # compute readsorption rate (U_SHE and pH independent)
    engsU = adjust_SHE_engs(mkin.engs0, U_SHE, mkin.echem)
    drt = make_rates(engsU, pH, [False]) 
    
    # evaluate flow (for one product hard-coded):
    #   i_diss > indicator which pr rates go to ketene-product
    #   p0 * rate_readsorption * empty sites
    flux0 = sum(pr[:i_diss]) - p0[0]*drt[i_rads]*ys[-1,0]
    # pr = [C_t, D_t]
    return(ts, ys, [flux0, *pr[i_diss:]])


def _normalization_flux(mkin, y0, U, pH, pCO, i_diss, i_rads, pr, conv):
    """
      helper-function to normalize production rate
      this function is left in to double check, a missing mass-conservation
      was in this model due to an error in the mkm-equations;
      
      could in principle to depreciated (but doesn't hurt neither)

    """
    # post-normalize flux created from extra source term + Warning
    # (a) flux with NO extra outer pressure
    ts00, ys00, pr00 = _eval_flux_mkm(mkin, y0, U, pH, [0.0, pCO], i_diss, i_rads)
    # (b) add Warning/message if outer flux > 0.01%
    rtot = sum(pr) / sum(pr00)
    if rtot > 1.0001 and rtot < 10.: # left in to double check
        logging.info("TRNSouter: Warning: increase in total rate through" + \
            "transport %.2e vs. %.2e (ratio %.3f)"%(sum(pr), sum(pr00), rtot-1.0))
    elif rtot > 10: # total rate increase > 10 counted as not converged
        logging.info("TRNSouter: Failed convergence: increase in total rate through" + \
            "transport %.2e vs. %.2e (ratio %.3f)"%(sum(pr), sum(pr00), rtot-1.0))
    # (c) renormalize pr (not ys -- although should)
    pr /= rtot
    # (d) give np.nan to non-converged values
    if conv > 0.5 or rtot > 10:
        pr[:] = np.nan
    return(pr, rtot)


def _warn_conv(conv, U, pH, Lx, diffmode, p0, c0):
    if conv > 1e-4:
        logging.info("TRNSouter: Warning:  %.2f V, pH=%.2f, Lx=%.2e, %s did not converge:\n"%(U,pH,Lx,diffmode) + \
            "conv = %.3e; p0=%.3e vs. c0(analytic)=%.3e"%(conv, p0[0], c0))


def compute_ketene_c0(flux, Lx, pH, roughness=1.0): #hardcoded ketene desorption
    """
      function for analytical surface concentration

    """
    # NOTE: conversions in cm2/s
    fluxM = convert_TOF2flux(flux, roughness)

    # note on conversions 1e4 (m2/s --> cm2/s); 1e3 (mol/m3 * s --> mol/l * s)
    c0a = get_analytical_conc0(fluxM, Lx, c_D_kt*1e4)
    c0b = get_analytical_conc0_diffrct(fluxM, Lx, c_D_kt*1e4, c_k_sol*1e3, pH)
    # output in ml/cm3 = 1e3 * std acitivity...
    c0a *= 1e3; c0b *= 1e3
    return({'diff':c0a, 'diff-rct':c0b})


def get_analytical_conc0(cgrd, Lx, D, cbulk=0):
    """
      helper-function for analytical surface concentration of ketene in
      diffusion-only case

    """
    ## just Fick's 1st law of diffuison: cgrd = D*dy/dx
    dy = (cgrd*Lx)/D
    y0 = dy + cbulk
    return(y0)


def get_analytical_conc0_diffrct(cgrd, Lx, D, kexp, ph, cbulk=0):
    """
      helper-function for analytical surface concentration of ketene in
      diffusion-reaction

    """
    ke = kexp * 10**(ph-14) # kexp in 1/(mol/l) pH in mol/l --> 1/s
    cph = np.sqrt(ke*D) # boundary == inf
    c0 = cgrd / cph
    return(c0)


