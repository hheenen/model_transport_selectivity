#!/usr/bin/env python

"""
This script executes the key function of `model_COR_Ac` to simulate the 
roughness, pH, and potential dependence of the acetate-selectivity model 
and evaluate and plot reproducing the plots in the manuscript and SI

HERE: Sensitivity of model

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import golden_ratio

# internal functions of model_COR_Ac to evaluate concentration profiles
from model_COR_Ac.simtools import potsweep_rgh_pH, make_facet_mkin_model_gpaw 
from model_COR_Ac.plotting import set_plotting_env, clrs, writefig


def sample_readsorption_energy():
    """
      helper-function to screen readsorption energy vs. selectivity

    """
    # predefined readsorption energy range
    Gads_range = np.array(np.arange(0.3, 0.81, 0.05).tolist()+[1.1])
    U_SHE = [-1.3]; pH_bulk = [14.0]; rghs = [15, 60]
    
    # prep eng-model
    mkin = make_facet_mkin_model_gpaw(100)
    br_ex = mkin.engs0[3][1] - mkin.engs0[3][0]

   
    # iterate energy range
    dat = {rgh:np.zeros((Gads_range.size, 2)) for rgh in rghs}
    for i in range(Gads_range.size):
        Gads = Gads_range[i]
        # adjust energies
        mkin.engs0[2][2] = -0.4 + br_ex - Gads
        mkin.engs0[3][2] = br_ex - Gads
        
        d = potsweep_rgh_pH(Us=U_SHE, pHs=pH_bulk, \
            rghs=rghs, Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', \
            mkin=mkin, savekey='eng_screen', tol_diff=1e-8)
        for rgh in rghs:
            dat[rgh][i,:] = (Gads, d[rgh][14.0]['SAcC2'][0])
    return(dat)


def plot_sel_vs_Gads(filename, dat, pts):
    """
      plot-function to plot selectivity vs. Gads
      
      Parameters
      ----------
      filename : str
        prefix of figure .pdf
      dat : dict
        simulated data
      pts : list
        DFT reference data

    """
    set_plotting_env(width=3.37,height=3.37/ golden_ratio,\
                   lrbt=[0.15,0.95,0.18,0.93],fsize=7.0)

    # colors - curves
    rclr = {15:'k', 60:'darkgray'}

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # plot selectivity curve for different roughnesses
    for rgh in dat:
        ax1.plot(dat[rgh][:,0], dat[rgh][:,1]*100, ls='-', color=rclr[rgh], label=r"$\rho$=%i"%rgh)

    # axis labels
    ax1.set_xlabel(r'$\Delta G^{\ddag}_{\mathrm{ads}}$ (eV)')
    ax1.set_ylabel(r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')

    # markers for DFT data - partly hard-coded
    dx = {110:[pts[110]-0.05, 85], 211:[pts[211]-0.05,85], 100:[pts[100]-0.05, 60], 111:[pts[111]-0.05, 40]}
    for surf in pts:
        ax1.axvline(pts[surf], color=clrs['orange'], ls='--')
        t = plt.text(*dx[surf], r"Cu(%i)"%surf, fontsize=6, color=clrs['orange'])
        t.set_bbox(dict(facecolor='white', edgecolor='white', pad=0.05))

    ax1.axvspan(0.475, 0.525, color=clrs['lightblue'], alpha=0.4)
    ax1.annotate('exp', xy=(0.5,10), xytext=(0.55, 20), size=6, color=clrs['deepblue'],\
        arrowprops={"arrowstyle":'-'})
   
    ax1.legend(loc=2,prop={'size':6},handlelength=0.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))

    writefig(filename ,folder='output') 


def sample_rgh_pH():
    """
      helper-function to screen readsorption energy vs. selectivity

    """
    mkin = make_facet_mkin_model_gpaw(100)
    
    U_SHE = [-1.3]

    # sample pH
    pH_range = np.around(np.arange(12.0, 14.6, 0.2),2)
    dpH = potsweep_rgh_pH(Us=U_SHE, pHs=pH_range, \
        rghs=[1], Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pH_rgh_screen', tol_diff=1e-8)
   
    # sample rgh
    rghs = np.array(list(range(1,10))+list(range(10,101,5)))
    drgh = potsweep_rgh_pH(Us=U_SHE, pHs=[14.0], \
        rghs=rghs, Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pH_rgh_screen', tol_diff=1e-8)

    return(dpH, drgh)


def plot_selectivity_vs_rgh_pH(filename, dph, drgh):
    """
      plot-function to plot selectivity vs. rgh and pH
      
      Parameters
      ----------
      filename : str
        prefix of figure .pdf
      dph : dict
        simulated data - pH
      drgh : dict
        simulated data - rgh

    """
    set_plotting_env(width=3.37,height=3.37/ golden_ratio,\
                   lrbt=[0.13,0.87,0.18,0.93],fsize=7.0)


    # colors
    rclr = {15:'k', 60:'darkgray'}

    fig, ax = plt.subplots(1,2) # first only pH

    # prepare data - pH
    dph_sel = np.array([[10**(pH-14), dph[1][pH]['SAcC2'][0]] for pH in dph[1]])
    dph_ckt = np.array([[10**(pH-14), dph[1][pH]['B_t'][0]] for pH in dph[1]])
    dph_cc2 = np.array([[10**(pH-14), dph[1][pH]['C_t'][0]] for pH in dph[1]])

    # plot data - pH
    ax[0].plot(dph_sel[:,0], dph_sel[:,1]*100., color=clrs['darkgray'])
    axa = ax[0].twinx()
    axa.plot(dph_ckt[:,0], dph_ckt[:,1], color=clrs['deepblue'], ls=':')
    axa.plot(dph_cc2[:,0], dph_cc2[:,1], color=clrs['darkred'], ls=':')
    ax[0].set_xscale('log')

    # prepare data - rgh
    drgh_sel = np.array([[rgh, drgh[rgh][14]['SAcC2'][0]] for rgh in drgh])
    drgh_ckt = np.array([[rgh, drgh[rgh][14]['B_t'][0]/rgh] for rgh in drgh])
    drgh_cc2 = np.array([[rgh, drgh[rgh][14]['C_t'][0]/rgh] for rgh in drgh])

    # plot data - pH
    ax[1].plot(drgh_sel[:,0], drgh_sel[:,1]*100., color=clrs['darkgray'])
    axb = ax[1].twinx()
    axb.plot(drgh_ckt[:,0], drgh_ckt[:,1], color=clrs['deepblue'], ls=':')
    axb.plot(drgh_cc2[:,0], drgh_cc2[:,1], color=clrs['darkred'], ls=':')

    # legend
    ax[1].plot(np.nan, np.nan, color=clrs['darkgray'], label=r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')
    ax[1].plot(np.nan, np.nan, color=clrs['deepblue'], ls=':', label=r'TOF (Ac$^-$)')
    ax[1].plot(np.nan, np.nan, color=clrs['darkred'], ls=':', label=r'TOF (C$_{2}$)')
    ax[1].legend(loc=7,prop={'size':6},handlelength=0.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))
    
    ax[0].set_ylim(0.,100.); ax[1].set_ylim(0.,100.)
    axa.set_ylim(0.0,0.06); axb.set_ylim(0.0,0.06)

    axa.yaxis.set_ticklabels([])
    ax[1].yaxis.set_ticklabels([])

    ax[0].set_ylabel(r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')
    axb.set_ylabel(r'TOF (1/(site $\cdot$ s))')
    ax[0].set_xlabel(r'[OH$^-$]')
    ax[1].set_xlabel(r'roughness $\rho$')

    plt.subplots_adjust(wspace=0.05)

    writefig(filename ,folder='output') 
    
    
def sample_polarization_curves_facets(include_ads=False):
    """
      helper-function to run model with exchanged energetics according to
      Cu(100), Cu(111), Cu(110), Cu(211) only critical energies are 
      maintained (see "make_facet_mkin_model_gpaw" in 
      "model_COR_Ac.mkm.model_surfrct"); optionally the adsorption energies
      are also adjusted Cu(100)

    """
    # range to sample for sensitivity analysis
    Us = np.arange(-1.7, -1.05, 0.1); pHs=[13.7]; rghs = [15]
    
    # preparation of energetics
    mkin100 = make_facet_mkin_model_gpaw(100); Eads100 = mkin100.engs0[0]
    dout = {}
    for surf in [111, 100, 110, 211]:
        mkin = make_facet_mkin_model_gpaw(surf)
        if include_ads:
            mkin.engs0[0] = Eads100
        
        dat = potsweep_rgh_pH(Us=Us, pHs=pHs, \
            rghs=rghs, Lx=[0.35e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='facet_screen2', tol_diff=1e-8)
        dout.update({surf:dat[15]})
    
    return(dout)


def plot_polarization_facet(filename, dat, adstag=False, legendtag=False):
    """
      plot-function to plot selectivity vs. potential
      for different facets-models
      
      Parameters
      ----------
      filename : str
        prefix of figure .pdf
      dat : dict
        simulated data
      adstag : bool
        if adding the tag about adsorption energies
      legendtag : bool
        if adding legend

    """
    set_plotting_env(width=3.37,height=3.37/ golden_ratio,\
                   lrbt=[0.15,0.95,0.15,0.98],fsize=7.0)

    # linestyles, colors
    lss = {111:':', 100:'-', 110:'--', 211:'-.'}
    pHs = list(dat[list(dat.keys())[0]].keys())
    pHclr = {13.7:clrs['deepblue'], 14.0:clrs['green'], 14.3: clrs['darkred']}

    fig, ax = plt.subplots() # first only pH
    for facet in dat:
        for pH in dat[facet]:
            d = dat[facet][pH]
            ax.plot(d['v'], np.array(d['SAcC2'])*100., ls=lss[facet], color=pHclr[pH])

    # legend
    for r in lss:
        ax.plot(np.nan, np.nan, ls=lss[r], color='k', label=r'Cu(%i)'%r)
    for pH in pHclr:
        if pH == 13.7:
            ax.plot(np.nan, np.nan, color=pHclr[pH], label=r'$\mathrm{pH_{bulk}}$=%.1f'%pH)
    if legendtag:
        ax.legend(loc=3,prop={'size':6},handlelength=1.5, frameon=False)

    if adstag:
        ax.annotate(r'Cu(100)-like $\Delta G\mathrm{_{ads}^{CO}}$', xy=(-1.35, 80))

    ax.set_ylim(0,100)

    ax.set_ylabel(r'$S\mathrm{_{Ac^-}^{C_{2}}}$ (%)')
    ax.set_xlabel(r'U vs. SHE (V)')

    writefig(filename ,folder='output') 


if __name__ == "__main__":

   #engs = {100:{'de':0.34, 'ad':0.71}, # SCCS
   #        111:{'de':0.05, 'ad':0.78},
   #        110:{'de':0.14, 'ad':0.43},
   #        211:{'de':0.22, 'ad':0.51}}
    engs = {100:{'de':0.36, 'ad':0.77}, # GPAW no water
            111:{'de':0.20, 'ad':1.05},
            110:{'de':0.37, 'ad':1.00},
            211:{'de':0.21, 'ad':0.52}}

    ##########################################
    ###  (1) sampling readsorption energy  ###
    dat = sample_readsorption_energy()
    plot_sel_vs_Gads("FigS11_readsorption_vs_Gads", dat, {k:engs[k]['ad'] for k in engs})
    ##########################################

    
    ##################################################
    ###  (2) sampling rghs and pH vs. selectivity  ###
    dph, drgh = sample_rgh_pH()
    plot_selectivity_vs_rgh_pH('FigS6_selectivity_vs_rgh_pH', dph, drgh)
    ##################################################


    ################################################################
    ###  (3) exchange facet thermodynamics - keep SDS energies  ####
    dat = sample_polarization_curves_facets()
    plot_polarization_facet('FigS13a_selectivity_facet', dat,\
        adstag=False, legendtag=False)
    dat = sample_polarization_curves_facets(include_ads=True)
    plot_polarization_facet('FigS13b_selectivity_facet', dat,
        adstag=True, legendtag=True)
    ################################################################


