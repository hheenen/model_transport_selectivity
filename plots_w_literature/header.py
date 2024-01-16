"""

plotting header for exp-plots

"""

import os, sys
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
sys.path.append("../.")
from sim_diffusion_selectivity import sample_data

sys.path.append("../experimental_reference")
from read_data import *
from model_plots import plot_xy, plot_xy_ax, plot_xy_ax_in, set_plotting_env, writefig, golden_ratio, clrs


def get_adsdes_engs(cd):
    dG = cd['dG']; ddG = cd['ddG']
    
    cH = cd['H']; kb = 8.6173e-5
    acE = kb*300*np.log(cH) # energy penalty to emulate activity coefficient
    if dG > 0:
      # dGdes = dG+ddG-acE; dGads = ddG
        dGdes = dG+ddG; dGads = max(ddG-acE, 0)
    else:
      # dGdes = ddG-acE; dGads = -1*dG+ddG
        dGdes = ddG; dGads = (-1*dG+ddG)-acE
    return(dGdes, dGads)


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
