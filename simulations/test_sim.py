#!/usr/bin/env python

import sys
sys.path.append('../.') # load path for tests

import numpy as np

from model_Sel_CORAc.simtools import potsweep_rgh_pH, make_facet_mkin_model_gpaw
from model_Sel_CORAc.tools import load_pickle_file


if __name__ == "__main__":
    mkin = make_facet_mkin_model_gpaw(100)
    Us = np.arange(-1.7, -1.05, 0.2); pHs=[14.3]; rghs = [15]
    Us = [-1.7, -1.4]; pHs=[14.3]; rghs = [15]

    dat = potsweep_rgh_pH(Us=Us, pHs=pHs, \
        rghs=rghs, Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='test_sim_gpaw', tol_diff=1e-8)
