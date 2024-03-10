#!/usr/bin/env python

import numpy as np
from model_Sel_CORAc.plotting import set_plotting_env, clrs, writefig
from scipy.constants import golden_ratio
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def plot_xy(filename, dat, ylabel, xlabel, ls_args, derr={}, \
    tag="", line=[], fsize=8.0, lsize=6.0, legpos=1, legcol=1, \
    csp=1.0, extra_leg=[], nolegend=False, title=''):
    """
      plotting function for xy-plots
      various parameters are taken as self-explanatory

    """
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05,\
               lrbt=[0.16,0.98,0.17,0.92],fsize=fsize)

    plot_xy_ax(dat, ylabel, xlabel, ls_args, derr,  \
        tag, line, lsize, legpos, legcol, csp, extra_leg, nolegend, title)

    writefig(filename ,folder='output', write_pdf=False, write_png=True) 


def plot_xy_ax(dat, ylabel, xlabel, ls_args, derr={},\
    tag="", line=[], lsize=6.0, legpos=1, legcol=1, \
    csp=1.0, extra_leg=[], nolegend=False, title=''):
    """
      helper function for plot_xy wrapping plt.figure environment 
      around plot_xy_ax_in

    """
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1 = plot_xy_ax_in(ax1, dat, ylabel, xlabel, ls_args, derr,\
        tag, line, lsize, legpos, legcol, csp, extra_leg, nolegend, title)

    return ax1


def plot_xy_ax_in(ax1, dat, ylabel, xlabel, ls_args, derr={},\
    tag="", line=[], lsize=6.0, legpos=1, legcol=1, \
    csp=1.0, extra_leg=[], nolegend=False, title=''):
    """
      helper function making actual plot

    """
    
    for k in dat:
        if k not in derr:
            derr.update({k:np.zeros(dat[k].shape)})

    for k in dat:
        xy = dat[k]
        err = derr[k]
        if err.sum() > 0:
            zz = 1
            if 'zorder' in ls_args[k]:
                zz = ls_args[k].pop('zorder')
            if len(extra_leg) == 0:
                ax1.errorbar(xy[:,0], xy[:,1], xerr=err[:,0], yerr=err[:,1], \
                    **ls_args[k], elinewidth = 0.5, capsize=1.0, zorder=zz, label=k)
            else:
                ax1.errorbar(xy[:,0], xy[:,1], xerr=err[:,0], yerr=err[:,1], \
                    **ls_args[k], elinewidth = 0.5, capsize=1.0, zorder=zz)
        elif len(extra_leg) == 0:
            ax1.plot(xy[:,0], xy[:,1], **ls_args[k], label=k)

        else:
            ax1.plot(xy[:,0], xy[:,1], **ls_args[k])

    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
   
    # make separat label dict
    if len(extra_leg) > 0:
        for kk in extra_leg:
            ax1.plot(np.nan, np.nan, **extra_leg[kk], label=kk)
    if not nolegend:
        ax1.legend(loc=legpos,ncol=legcol, columnspacing=csp, \
            prop={'size':lsize},handlelength=1.5, frameon=False)
   
    if len(line) > 0:
        if type(line) == list:
            for l in line:
                ax1.plot(l[:,0], l[:,1], color='r', ls='--', lw=0.5)
        else:
            ax1.plot(line[:,0], line[:,1], color='r', ls='--', lw=0.5)

    ax1.annotate(tag, xy=(0.25, 0.9), xycoords='axes fraction', size=8)
    ax1.annotate(title, xy=(0.45, 1.03), xycoords='axes fraction', size=10, \
        horizontalalignment='center')

    return ax1

