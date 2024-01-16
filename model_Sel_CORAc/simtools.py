"""
This module contains functions which serve as useful wrappers to 
conduct simulations with the multiscale module

"""

import sys, os
import numpy as np
import logging
from time import time
from datetime import datetime

from model_Sel_CORAc.tools import load_pickle_file, write_pickle_file, _make_hash
from model_Sel_CORAc.transport import solve_MS_model

# also load wrapper for simulations - to import 
from model_Sel_CORAc.mkm.model_surfrct import make_facet_mkin_model_gpaw


def potsweep_rgh_pH(Us, pHs, rghs, Lx, diffmode, mkin, savekey='stddat', tol_diff=None):
    """
      wrapper to simulate steady-state current densities and from this
      selectivities at a range for pH and potential
      
      Parameters
      ----------
      Us : list/array
        values of applied potential vs. SHE
      pHs : list/array
        values of applied pH
      rghs : list/array
        values of roughness of catalyst
      Lx : float
        diffusion length in cm 
      diffmode : str
        keyword for mode of ketene transport (diff / diff-rct / diff-rct-num)
      mkin : mkm-object
        pre-created mkm object (see model_Sel_CORAc.mkm.model_surfrct)
      savekey : str
        name of files to save results of simulation
      tol_diff : float
        tolerance for numerical solver to converge

      Returns
      -------
      dout : dict
        dictionary sorted after rgh: pH: with pot ('v'), 
        current density Ac ('B_t'), other C2 ('C_t'), Ac selectivity ('SAcC2'),
        concentration profiles for solution
    
    """
    # for saving data - save output in same folder
    pklfile = 'dat_%s.pkl'%savekey
    
    # use root logger to track convergence of routine
    logfile = 'dat_%s.log'%savekey
    
    # use info to be printed to stdout | base is DEBUG
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(level=logging.DEBUG,
                    handlers=[logging.FileHandler(logfile), 
                              stream_handler])
    
    logging.info(">>> run potsweep_rgh_pH on %s"%(\
        datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
    
    # identifier in pkl-file: clear paring of same mkm and transport parameters
    pkey = _make_prop_key(Lx, diffmode, mkin)
    
    # load pkl-file if exists -- update dout with pkey-matching data
    dsv = {pkey:{}}; dout = dsv[pkey]
    if os.path.isfile(pklfile):
        dsv = load_pickle_file(pklfile)
        if pkey not in dsv:
            dsv.update({pkey:{}})
        dout = dsv[pkey]

    # iterate through rgh, pH and U
    for rgh in rghs:
        if rgh not in dout:
            dout.update({rgh:{}})
        for pH in pHs:
            pH = np.around(pH,1)
            if pH not in dout[rgh]:
                dout[rgh].update({pH:{'v':[], 'B_t':[], 'C_t':[], 'SAcC2':[], 'Cnum':[], 'p_norm':[], 'p_raw':[]}})
            dat = dout[rgh][pH]
            for U in Us:
                U = np.around(U,2)
                if U not in dout[rgh][pH]['v']:
                    logging.info(">>> Beginning for rgh %.0f, pH %.1f, %.2f V_SHE"%(rgh, pH, U))
                    t0 = time()
                    ts, ys, pr, p0, p0_raw, Cx = solve_MS_model(mkin, U, pH, Lx, diffmode, i_diss=2, i_rads=7, roughness=rgh, tol_diff=tol_diff)
                    logging.info('converging automatically in %.1f s'%(time()-t0))
                    # store data
                    dat['v'].append(U); dat['B_t'].append(pr[0]); dat['C_t'].append(pr[1])
                    dat['SAcC2'].append(pr[0]/sum(pr))
                    dat['Cnum'].append(Cx)
                    dat['p_norm'].append(p0) # p_norm and p_raw (should always be the same) for debugging
                    dat['p_raw'].append(p0_raw) # pressures (=activities) are surface concentrations from mkm
                    dout[rgh][pH] = _sort_pot(dat)
                    write_pickle_file(pklfile, dsv)
    dout = {rgh:{pH:dout[rgh][pH] for pH in pHs} for rgh in rghs}
    return(dout)


def _sort_pot(dat):
    """
      helper-function to sort potentials -- saved in ascending order

    """
    Us = dat['v']
    ind = np.argsort(Us)
    for k in dat:
        dat[k] = [dat[k][i] for i in ind]
    return(dat)


def _make_prop_key(Lx, diffmode, mkin):
    """
      helper-function to create a key-identifier

    """
    # make key via to identify properties of simulation
    enghsh = _make_hash(np.array(mkin.engs0))
    key = '%.3e_%.3e_%s_%s'%(*Lx,diffmode,enghsh)
    return(key)
            

