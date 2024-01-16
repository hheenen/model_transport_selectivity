#!/usr/bin/env python

"""
This script executes the key function of `model_COR_Ac` to simulate the 
roughness, pH, and potential dependence of the acetate-selectivity model 
and evaluate and plot reproducing the plots in the manuscript and SI

HERE: main script to simulate current densities and from there selectivity

Note: the numerical solution of the solver is not very efficient, the 
script may run up to 20 minutes for all solutions on a single core 

"""

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from model_COR_Ac.simtools import potsweep_rgh_pH, make_facet_mkin_model_gpaw 
from model_COR_Ac.plotting import set_plotting_env, clrs, writefig
from model_COR_Ac.transport.flux_conversion import f2j



def sample_polarization_curves():
    """
      helper-function to sweep through potentials, pH and roughness via
      potsweep_rgh_pH

    """
    mkin = make_facet_mkin_model_gpaw(100)
    rghs = [10, 13, 65] #5, 15, 65]
    pHs=[13.7, 14.0, 14.3]
    Us = np.arange(-1.7, -1.05, 0.1)


    dat = potsweep_rgh_pH(Us=Us, pHs=pHs, \
        rghs=rghs, Lx=[0.15e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pol_screen_gpaw', tol_diff=1e-8)
    
    # rescale according to loading NS of 0.5 mg/cm2 (roughness=10)
    for pH in dat[10]:
        dat[10][pH]['B_t'] = (np.array(dat[10][pH]['B_t'])* 0.5).tolist()
        dat[10][pH]['C_t'] = (np.array(dat[10][pH]['C_t'])* 0.5).tolist()

    return(dat)


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
    ax[0].set_ylim(0.02, 1261.893951145012)
    #[x.set_xlim(-1.73, -1.07) for x in ax]

    ax[0].set_ylabel(r'$j\mathrm{_{geo}^{C_{2}}}$ (mA/cm$^2$)')
    ax[0].yaxis.set_label_coords(-0.14,0.5)
    ax[1].set_ylabel(r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')
    ax[1].set_xlabel(r'U vs. SHE (V)')

 #  ax[0].annotate(r'(a)', xy=(0.01, 0.96), xycoords='figure fraction', size=11)
    ax[0].annotate(r'Theory', xy=(0.50, 0.965), xycoords='figure fraction', size=9)

    ax[1].annotate(r'pH', xy=(0.04, 0.20), xycoords='axes fraction', size=7)
    ax[1].annotate(r'', xy=(0.10, 0.39), xytext=(0.10,0.08), xycoords='axes fraction', size=6, 
       arrowprops={"arrowstyle":'->', "color":'k', 'lw':1.0})
    ax[1].annotate(r'$\rho$', xy=(0.145, 0.20), xycoords='axes fraction', size=7)
    ax[1].annotate(r'', xy=(0.18, 0.38), xytext=(0.18,0.07), xycoords='axes fraction', size=6, 
       arrowprops={"arrowstyle":'<-', "color":'k', 'lw':1.0})

    plt.subplots_adjust(hspace=0.05)

    writefig(filename ,folder='output', write_png=True) 


def _prep_pol_dat(fdat):
    """
      helper-function to transform data flux to current densities and 
      selectivity

    """
    # rename keys
    jd = {}; fed = {}
    for key in fdat:
        d = {'Ac':np.array([fdat[key]['v'], f2j(fdat[key]['B_t'],4)]).T, \
            'C2+':np.array([fdat[key]['v'], f2j(fdat[key]['C_t'],8)]).T}
        jd.update({key:d})
        fe = {'Ac':np.array([fdat[key]['v'], (f2j(fdat[key]['B_t'],4)/(f2j(fdat[key]['B_t'],4)+f2j(fdat[key]['C_t'],8))*100.)]).T, \
             'C2+':np.array([fdat[key]['v'], (f2j(fdat[key]['C_t'],8)/(f2j(fdat[key]['B_t'],4)+f2j(fdat[key]['C_t'],8))*100.)]).T}
        fed.update({key:fe})
    return(jd, fed)



if __name__ == "__main__":
    
    ################################################
    #### sampling pot and rghs and pH - Fig. 4a ####
    dat = sample_polarization_curves()
    plot_polarization_sel('Fig4a_selectivity_pol', dat)
    ################################################
   

