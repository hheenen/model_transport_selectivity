"""

plotting header for plots

"""

# general imports
import os, sys
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

# package imports
sys.path.append("../.")
from model_Sel_DRAR.sim_diffusion_selectivity import sample_data
from model_Sel_DRAR.model_plots import plot_xy, plot_xy_ax, set_plotting_env, \
                                        writefig, golden_ratio, clrs

# imports experimental data
sys.path.append("../experimental_reference")
from read_data import *


def get_adsdes_engs(cd):
    """
      helper-function to convert barriers to thermodynamics as barriers 
      (see SI)

    """
    dG = cd['dG'][0]; ddG = cd['ddG_des'][0]
    
    cH = cd['cH'][0]; kb = 8.6173e-5
    acE = kb*300*np.log(cH) # energy penalty to emulate activity coefficient
    if dG > 0:
      # dGdes = dG+ddG-acE; dGads = ddG
        dGdes = dG+ddG; dGads = max(ddG-acE, 0)
    else:
      # dGdes = ddG-acE; dGads = -1*dG+ddG
        dGdes = ddG; dGads = (-1*dG+ddG)-acE
    return(dGdes, dGads)


def run_model_for_example(cdict, sdict, okey='roughness'):
    """
      helper-function to make full plot
      note: okey could be smartly seeing what is a range in sdict...

    """
    out_plot = {}
    for k in sdict:
        # prep energies for model input
        dGdes , dGads = get_adsdes_engs(cdict[k])

        datafile = "model_data_examples_%s.pkl"%k
        # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
        dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=cdict[k]['ddG_red'], \
            Ds=cdict[k]['cD'], Lxs=[50e-4], **sdict[k])
        # save roughness vs. index 
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        iout = ['i_model', 'eng_des', 'eng_ads', 'eng_red', \
            'U_SHE', 'Dx', 'Lx', 'roughness']
        out_plot.update({k:np.array([dat[:,iout.index(okey)], sel]).T})

    return out_plot


def add_sketch(ax, a1, a2, a3, dxy=(0.0,0.0)):
    dx, dy = dxy
    ax.annotate(a1, xy=(0.46+dx,0.5+dy), xycoords='axes fraction', ha='right', size=9)
    ax.annotate(a2, xy=(0.65+dx,0.5+dy), xycoords='axes fraction', ha='left', size=9)
    ax.annotate(a3, xy=(0.52+dx,0.71+dy), xycoords='axes fraction', ha='left', size=9)
    ax.annotate(r'', xy=(0.65+dx,0.53+dy), xytext=(0.47+dx,0.53+dy), xycoords='axes fraction',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    ax.annotate(r'', xy=(0.58+dx,0.7+dy), xytext=(0.50+dx,0.53+dy), xycoords='axes fraction',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.5"))
    
def add_sketch2(ax, a1, a2, a3, dxy=(0.0,0.0)):
    dx, dy = dxy
    ax.annotate(a1, xy=(0.05+dx,0.88+dy), xycoords='axes fraction', ha='left', size=9)
    ax.annotate(a2, xy=(0.05+dx,0.6+dy), xycoords='axes fraction', ha='left', size=9)
    ax.annotate(a3, xy=(0.15+dx,0.71+dy), xycoords='axes fraction', ha='left', size=9)
    ax.annotate(r'', xy=(0.1+dx,0.65+dy), xytext=(0.1+dx,0.86+dy), xycoords='axes fraction',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    ax.annotate(r'', xy=(0.18+dx,0.74+dy), xytext=(0.1+dx,0.84+dy), xycoords='axes fraction',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.5"))

