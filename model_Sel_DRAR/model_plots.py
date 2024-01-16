#!/usr/bin/env python

import numpy as np
from model_Sel_CORAc.plotting import set_plotting_env, clrs, writefig
from scipy.constants import golden_ratio
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def plot_polarization_sel(filename, dat):
    """
      plot-function to plot sweep through potentials, pH and roughness,
      two plots ax[0] --> partial current density
                ax[1] --> selectivity

    """
    set_plotting_env(width=3.37,height=3.37,\
                   lrbt=[0.17,0.98,0.11,0.95],fsize=9.0)

    # linestyles, colors and markers
    lss = {10:'-', 13:':', 65:'--'}
    #lss = {5:'-', 15:':', 65:'--'}
    pHs = list(dat[list(dat.keys())[0]].keys())
    pHclr = {13.7:clrs['deepblue'], 14.0:clrs['green'], 14.3: clrs['darkred']}
    
    # prepare data - pH
    jd = {r:_prep_pol_dat(dat[r])[0] for r in dat}
    
    fig, ax = plt.subplots(2,1) # first only pH

    # arrows and labels in plot -- hard-coded
    prps = {'color':clrs['orange'], 'size':7}
    aprps = {"arrowstyle":'->', 'color':clrs['orange']}
    ax[1].annotate("kinetic limited", xy=(-1.33,84), rotation=12, **prps)
    ax[1].annotate("(SDS-surf)", xy=(-1.29,72), rotation=13, **prps)
    ax[1].annotate("", xy=(-1.35,76), xytext=(-1.13,92), arrowprops=aprps)
    ax[1].annotate("transport limited", xy=(-1.60,82.5), rotation=-10., **prps)
    ax[1].annotate("(SDS-sol)", xy=(-1.56,71), rotation=-11., **prps)
    ax[1].annotate("", xy=(-1.63,89), xytext=(-1.38,75), arrowprops=aprps)
    
    # iterate roughness & pH in data
    for rgh in dat:
        for pH in dat[rgh]:
            d = dat[rgh][pH]
            ja = jd[rgh][pH]['Ac']; jc = jd[rgh][pH]['C2+']; jt = deepcopy(ja); jt[:,1] += jc[:,1]
            ax[0].plot(jt[:,0], jt[:,1], ls=lss[rgh], color='k') # current
            ax[1].plot(d['v'], np.array(d['SAcC2'])*100., ls=lss[rgh], color=pHclr[pH]) # selectivity

    # legend
    for r in lss:
        ax[0].plot(np.nan, np.nan, ls=lss[r], color='k', label=r'$\rho$=%i'%r)
    for pH in pHclr:
        if pH > 13.0:
            ax[0].plot(np.nan, np.nan, color=pHclr[pH], label=r'$\mathrm{pH_{bulk}}$=%.1f'%pH)
    ax[0].legend(loc=3,prop={'size':6},handlelength=1.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))

    # format axes
    ax[0].xaxis.set_ticklabels([])
    ax[1].set_ylim(0,100)
    ax[0].set_yscale('log')

    ax[0].set_ylabel(r'$j\mathrm{_{geo}^{C_{2}}}$ (mA/cm$^2$)')
    ax[0].yaxis.set_label_coords(-0.14,0.5)
    ax[1].set_ylabel(r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')
    ax[1].set_xlabel(r'U vs. SHE (V)')

    ax[0].annotate(r'(a)', xy=(0.01, 0.96), xycoords='figure fraction', size=11)
    ax[0].annotate(r'Theory', xy=(0.50, 0.965), xycoords='figure fraction', size=9)

    ax[1].annotate(r'pH', xy=(0.04, 0.20), xycoords='axes fraction', size=7)
    ax[1].annotate(r'', xy=(0.10, 0.39), xytext=(0.10,0.08), xycoords='axes fraction', size=6, 
       arrowprops={"arrowstyle":'->', "color":'k', 'lw':1.0})
    ax[1].annotate(r'$\rho$', xy=(0.145, 0.20), xycoords='axes fraction', size=7)
    ax[1].annotate(r'', xy=(0.18, 0.38), xytext=(0.18,0.07), xycoords='axes fraction', size=6, 
       arrowprops={"arrowstyle":'<-', "color":'k', 'lw':1.0})

    plt.subplots_adjust(hspace=0.05)

    writefig(filename ,folder='output', write_png=True) 


def plot_selectivity_maps(filename, v, dat, xlabel, ylabel):
    z = np.array([]); xy = np.array([])
    for pH in dat:
        # xy
        y = np.ones(v.size)*pH
        xy0 = np.vstack((v,y)).T
        if xy.size == 0:
            xy = xy0
        else:
            xy = np.vstack((xy, xy0))
        z = np.hstack((z, dat[pH]))
    
    _plot_map(filename, xy, z, \
        ylabel=ylabel, xlabel=xlabel, zlabel=r'Selectivity $r_{\mathrm{EC}}$',\
        showmarkers=False)


def _plot_map(filename, xy, z, ylabel, xlabel, zlabel, showmarkers=True, tag="", line=[]):
    set_plotting_env(width=3.37,height=3.37/ golden_ratio,\
               lrbt=[0.16,0.9,0.16,0.95],fsize=8.0)
    
    # first make a scatter with continuous color map
    x = np.unique(xy[:,0]); y = np.unique(xy[:,1])
    zxy = np.zeros((y.size,x.size))
    for i in range(y.size):
        for j in range(x.size):
            indy = np.where(xy[:,1] == y[i])[0]
            ind = indy[np.where(xy[indy,0] == x[j])[0]]
            zxy[i,j] = z[ind]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    CS = ax1.contourf(x,y,zxy,8,vmin=z.min(),vmax=z.max(),cmap=plt.cm.GnBu)#jet)
    cbr = plt.colorbar(CS)
    
    # add points with value color to check interpolation
    cm = plt.get_cmap('GnBu')
    cNorm = colors.Normalize(vmin=z.min(),vmax=z.max())
    clrmap = plt.cm.ScalarMappable(norm=cNorm,cmap=cm)
    if showmarkers:
        ax1.scatter(xy[:,0],xy[:,1],marker='x',s=0.8,color='k')
        for i in range(xy[:,0].size):
            ax1.plot(xy[i,0], xy[i,1], marker='o', markersize=5.2, color=clrmap.to_rgba(z[i]))
    
    ax1.set_xlim(x.min(),x.max())
    ax1.set_ylim(y.min(),y.max())
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    cbr.ax.get_yaxis().labelpad = 12
    cbr.ax.set_ylabel(zlabel, rotation=270.)
   
    if len(line) > 0:
        if type(line) == list:
            for l in line:
                ax1.plot(l[:,0], l[:,1], color='r', ls='--', lw=0.5)
        else:
            ax1.plot(line[:,0], line[:,1], color='r', ls='--', lw=0.5)

    ax1.annotate(tag, xy=(0.2, 0.2), xycoords='axes fraction', size=6)

    writefig(filename ,folder='output', write_pdf=False, write_png=True) 
    

def plot_xy(filename, dat, ylabel, xlabel, ls_args, derr={},  tag="", line=[], fsize=8.0, lsize=6.0, legpos=1, legcol=1, csp=1.0, extra_leg=[], nolegend=False, title=''):
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05,\
               lrbt=[0.16,0.98,0.17,0.92],fsize=fsize)

    plot_xy_ax(dat, ylabel, xlabel, ls_args, derr,  tag, line, fsize, lsize, legpos, legcol, csp, extra_leg, nolegend, title)

    writefig(filename ,folder='output', write_pdf=False, write_png=True) 


def plot_xy_ax(dat, ylabel, xlabel, ls_args, derr={},  tag="", line=[], fsize=8.0, lsize=6.0, legpos=1, legcol=1, csp=1.0, extra_leg=[], nolegend=False, title=''):
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1 = plot_xy_ax_in(ax1, dat, ylabel, xlabel, ls_args, derr,  tag, line, fsize, lsize, legpos, legcol, csp, extra_leg, nolegend, title)

    return ax1


def plot_xy_ax_in(ax1, dat, ylabel, xlabel, ls_args, derr={},  tag="", line=[], fsize=8.0, lsize=6.0, legpos=1, legcol=1, csp=1.0, extra_leg=[], nolegend=False, title=''):
    
    for k in dat:
        if k not in derr:
            derr.update({k:np.zeros(dat[k].shape)})

    for k in dat:
        xy = dat[k]
        err = derr[k]
        if err.sum() > 0:
            ax1.errorbar(xy[:,0], xy[:,1], xerr=err[:,0], yerr=err[:,1], **ls_args[k], label=k, elinewidth = 0.8, capsize=2.0)
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
        ax1.legend(loc=legpos,ncol=legcol, columnspacing=csp ,prop={'size':lsize},handlelength=1.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))
   
    if len(line) > 0:
        if type(line) == list:
            for l in line:
                ax1.plot(l[:,0], l[:,1], color='r', ls='--', lw=0.5)
        else:
            ax1.plot(line[:,0], line[:,1], color='r', ls='--', lw=0.5)

    ax1.annotate(tag, xy=(0.25, 0.9), xycoords='axes fraction', size=8)
    ax1.annotate(title, xy=(0.45, 1.03), xycoords='axes fraction', size=10, horizontalalignment='center')

    return ax1
