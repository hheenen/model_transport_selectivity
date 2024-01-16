#!/usr/bin/env python

"""
This script executes the key function of `model_COR_Ac` to simulate the 
roughness, pH, and potential dependence of the acetate-selectivity model 
and evaluate and plot reproducing the plots in the manuscript and SI

HERE: Surface concentration behavior and profiles

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.constants import golden_ratio

# import wrapper from other simulation
from sim_current_selectivity import sample_polarization_curves, _prep_pol_dat

# internal functions of model_COR_Ac to evaluate concentration profiles
from model_COR_Ac.simtools import make_facet_mkin_model_gpaw 
from model_COR_Ac.transport import solve_MS_CO_analytical, \
                        estimate_convective_pH_analytical
from model_COR_Ac.transport.flux_conversion import f2j, convert_TOF2flux
from model_COR_Ac.plotting import set_plotting_env, clrs, writefig
from model_COR_Ac.tools import load_pickle_file, write_pickle_file


def plot_conc_profile_ketene(filename, x, c, pot, jac=[]):
    """
      plot-function to plot concentration profile of ketene
      
      Parameters
      ----------
      filename : str
        prefix of figure .pdf
      x : list/array
        x-axis in cm
      c : list of arrays
        concentration profiles
      pot : list
        list of potentials
      jac : array
        array of corresponding potential - Ac current pairs

    """
    set_plotting_env(width=3.37,height=3.37 / golden_ratio,\
                   lrbt=[0.17,0.98,0.17,0.98],fsize=9.0)

    # colors
    colors = plt.cm.jet(np.linspace(0,1,len(pot)))
    vclr = {pot[i]:colors[i] for i in range(len(pot))}
    # overwrite colors
    vclr = {-1.2:[1., 0.276, 0.], -1.4:[0., 0.486, 0.188], -1.6:[0., 0.6, 1.]}

    # plot concentration profiles
    fig, ax = plt.subplots()
    ipot = [-1.6, -1.4, -1.2] # hardcoded sorter
    dx = {-1.2:-0.0, -1.4:-0.14, -1.6:0.0}
    Is = []
    for u in ipot:
        i = pot.index(u)
        ax.plot(x*1e4, c[i][0,:]*1e-3*1e3, lw=2.0, color=vclr[pot[i]])
        
        kp = np.cumsum(c[i][0,:])/c[i][0,:].sum()
        iv = np.where(kp > 0.9)[0][0]; xv = x[iv]*1e4; yv = c[i][0,iv]
        #ax.axvline(xv, ls=':', color=vclr[pot[i]])
        ax.plot([xv, xv], [yv-0.002, yv+0.008], ls=(0,(1,1)), color=vclr[pot[i]])
        ax.annotate(r'10%', xy=(xv, yv+0.006), xytext=(xv+0.04+dx[u], yv+0.006), size=6, \
            color=vclr[pot[i]], arrowprops={"arrowstyle":'-', "color":vclr[pot[i]]})
        Is.append(i) # for fill_between

    for i in range(len(Is)-1):
        ax.fill_between(x*1e4, c[Is[i]][0,:], c[Is[i+1]][0,:], color=vclr[pot[Is[i]]], alpha=0.3, zorder=0)
    
    # format axes + legend
    ax.set_xlabel(r'distance from catalyst $x$ ($\mu m$)')
    ax.set_ylabel(r'[H$_2$CCO] (mmol/l)')

    ax.set_ylim(ax.get_ylim()[0], 0.055)
    ax.set_xlim(-0.02, 0.82)
    
    lstr = {v:r'%.1f $U_{\mathrm{SHE}}$'%v for v in ipot}
    if len(jac) > 0:
        lval = {v:jac[np.where(jac[:,0] == v)[0], 1][0] for v in ipot}
        print(lval)
        #prec = {v:3-len("%i"%lval[v]) for v in lval}
        prec = {-1.2:2, -1.4:1, -1.6:0} # hardcoded for paper -- odd looking decimals otherwise 
        istr = {v:"{:.{prec}f}".format(lval[v], prec=prec[v]) for v in lval}
        lstr = {v:'$j_{\mathrm{Ac^-}}$ = %s mA cm$^{-2}$'%(istr[v]) for v in lval}
    for v in lval:
        ax.plot(np.nan, np.nan, color=vclr[v], label=lstr[v])
    ax.legend(loc=5,prop={'size':7},handlelength=0.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))
    
    ax.annotate(r'(a)', xy=(0.0, 0.93), xycoords='figure fraction', size=11)
    
    # cartoon
    ax.annotate(r'', xy=(0.12, 0.025), xytext=(0.08,0.052),
            arrowprops={"arrowstyle":'->', "color":'darkgray'})
    ax.annotate(r'H$_2$CCO + OH$^- \rightarrow$ Ac$^-$', xy=(0.12,0.045), color='darkgray', size=8, weight='bold')


    writefig(filename ,folder='output', write_png=True) 


def sample_local_pH_convection():
    """
      helper-function to compute the local pH as obtained from our model

    """
    # local pH vs. current
    dat = {}
    pkl_file = 'dat_local_pH_convection.pkl'
    if os.path.isfile(pkl_file):
        dat = load_pickle_file(pkl_file)

    mkin = make_facet_mkin_model_gpaw(100)
    Lx = 0.35e-4 # diffusion length CO
    roughness = 15 # roughness doesn't really matter

    # solve coupled kinetic and CO/convection model
    # for different pH and potentials
    for pH in [7.0, 13.7, 14.0, 14.3]:
        if pH not in dat:
            dat.update({pH:[]})
            for U in np.arange(-1.8,-1.1,0.05):
                # (1) solve reactant
                ts, ys, pr, pCO = solve_MS_CO_analytical(mkin, U, pH, Lx, roughness)
                # (2) solve local pH
                pH_local2 = estimate_convective_pH_analytical(pr, pH, roughness)
                f_OH = convert_TOF2flux(pr[0]*4 + sum(pr[1:])*8, roughness)
                
                j = f2j(pr[0], 4)*roughness + f2j(pr[1], 8)*roughness
                dat[pH].append([U, j, pH, pH_local2])
            dat[pH] = np.array(dat[pH])
            write_pickle_file(pkl_file, dat)
    return(dat)


def plot_local_pH_convection(filename, dat):
    """
      plot-function to plot local pH for different conditions

    """
    set_plotting_env(width=3.37,height=3.37 / golden_ratio,\
                   lrbt=[0.15,0.95,0.16,0.98],fsize=7.0)
   
    # colors
    pHs = list(dat.keys())
    pHclr = {7.0:clrs['darkgray'], 13.0:clrs['azurblue'], 
            13.7:clrs['deepblue'], 14.0:clrs['green'], 14.3: clrs['darkred']}

    fig, ax = plt.subplots()
    for pH in dat:
        # local pH
        ax.plot(dat[pH][:,1], dat[pH][:,3], color=pHclr[pH])
        # bulk pH
        ax.plot([dat[pH][0,1], dat[pH][-1,1]], [pH,pH], color=pHclr[pH], ls=':', lw=1)
    ax.plot([0.0, dat[7][-1,1]], [7, dat[7][-1,3]], color=pHclr[7])

    # legend + axes
    for pH in dat:
        ax.plot(np.nan, np.nan, color=pHclr[pH], label=r'$\mathrm{pH_{bulk}}$=%.1f'%pH)
    ax.legend(loc=4,prop={'size':6},handlelength=0.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))
    
    ax.set_ylabel(r'local pH')
    ax.set_xlabel(r'$j\mathrm{_{geo}}$ (mA/cm$^2$)')

    writefig(filename ,folder='output') 


def get_surf_conc():
    """
      helper-function to obtain surface concentration from sampled data

    """
    # get data for ketene and OH (CO is not saved)
    dat = sample_polarization_curves()
    
    # ketene i = 0, OH i = 1
    csurf_kt = {rgh:{} for rgh in dat}
    csurf_oh = {rgh:{} for rgh in dat}
    for rgh in dat:
        for pH in dat[rgh]:
            v = dat[rgh][pH]['v']
            # only take first occurrence (on axis=1)
            ckt = [dat[rgh][pH]['Cnum'][ii][0,0] for ii in range(len(dat[rgh][pH]['Cnum']))]
            coh = [dat[rgh][pH]['Cnum'][ii][1,0] for ii in range(len(dat[rgh][pH]['Cnum']))]
            csurf_kt[rgh].update({pH:np.array([v,ckt]).T})
            csurf_oh[rgh].update({pH:np.array([v,coh]).T})
    
    # compute CO (not saved in general routine)
    pklfile = 'dat_pCO_surf.pkl'
    if os.path.isfile(pklfile):
        csurf_co = load_pickle_file(pklfile)
    else:
        mkin = make_facet_mkin_model_gpaw(100)
        csurf_co = {rgh:{pH:[] for pH in dat[rgh]} for rgh in dat}
        for rgh in dat:
            for pH in dat[rgh]:
                v = dat[rgh][pH]['v']
                for U in v:
                    ts, ys, pr, pCO = solve_MS_CO_analytical(mkin, U, pH, 0.35e-4, rgh)
                    csurf_co[rgh][pH].append([U, pCO])
                csurf_co[rgh][pH] = np.array(csurf_co[rgh][pH])
        write_pickle_file(pklfile, csurf_co)
    return(csurf_kt, csurf_oh, csurf_co)


def plot_surf_concentrations(filename, cco, coh, ckt):
    """
      plot-function to plot surface concentrations of 
      CO (cco input, ax[0]) OH (coh input, ax[1]) and 
      ketene (ckt input, ax[2]) for different conditions

    """
    set_plotting_env(width=3.37*2,height=3.37/ golden_ratio,\
                   lrbt=[0.07,0.99,0.17,0.98],fsize=9.0)

    # linestyles and colors
    lss = {10:'-', 13:':', 65:'--'}
    pHclr = {13.7:clrs['deepblue'], 14.0:clrs['green'], 14.3: clrs['darkred']}
    
    # plot sampled surface concentrations
    fig, axes = plt.subplots(1,3)
    for rgh in cco:
        for pH in cco[rgh]:
            axes[0].plot(cco[rgh][pH][:,0], cco[rgh][pH][:,1], ls=lss[rgh], color=pHclr[pH])
            axes[1].plot(coh[rgh][pH][:,0], coh[rgh][pH][:,1], ls=lss[rgh], color=pHclr[pH])
            axes[2].plot(ckt[rgh][pH][:,0], ckt[rgh][pH][:,1], ls=lss[rgh], color=pHclr[pH])


    # annotations
    axes[0].annotate(r'$c\mathrm{^{surf}_{CO}}$', xy=(0.7, 0.25), xycoords='axes fraction', size=11)
    axes[1].annotate(r'$c\mathrm{^{surf}_{OH^-}}$', xy=(0.7, 0.25), xycoords='axes fraction', size=11)
    axes[2].annotate(r'$c\mathrm{^{surf}_{H_2CCO}}$', xy=(0.7, 0.25), xycoords='axes fraction', size=11)
    
    axes[0].annotate(r'(a)', xy=(0.002, 0.94), xycoords='figure fraction', size=11)
    axes[1].annotate(r'(b)', xy=(0.35, 0.94), xycoords='figure fraction', size=11)
    axes[2].annotate(r'(c)', xy=(0.68, 0.94), xycoords='figure fraction', size=11)

    # axes formating 
    axes[0].set_ylabel(r'$c^{\mathrm{surf}}$ (mmol/l)')
    [axes[i].set_xlabel(r'U vs. SHE (V)') for i in range(len(axes))]
    plt.subplots_adjust(wspace=0.26)
    [axes[i].set_ylim(0.0, axes[i].get_ylim()[1]) for i in range(len(axes))]
    
    # legend
    for r in lss:
        axes[2].plot(np.nan, np.nan, ls=lss[r], color='k', label=r'$\rho$=%i'%r)
    for pH in pHclr:
        if pH > 13.0:
            axes[2].plot(np.nan, np.nan, color=pHclr[pH], label=r'$\mathrm{pH_{bulk}}$=%.1f'%pH)
    axes[2].legend(loc=1,prop={'size':6},handlelength=1.5, frameon=False)#, bbox_to_anchor=(1.0, 0.1))

    axes[2].yaxis.set_major_locator(MultipleLocator(0.02))

    writefig(filename ,folder='output', write_png=True) 


if __name__ == "__main__":

    ######################################
    ####       concentrations         ####
    ######################################

    #######################################################
    ### (1) H2CCO concentration profile - rgh=15, pH=14 ###
    d = sample_polarization_curves()
    c = d[13][14]['Cnum'] # concentration data
    jd, fed = _prep_pol_dat(d[13])
    j_ac = jd[14]['Ac']
    x = np.linspace(0.0,1e-4,c[0][0,:].size) # x in Lx
    plot_conc_profile_ketene('Fig5a_ketene_profile', x, c, d[13][14]['v'], j_ac)
    #######################################################
    quit()

    ##############################################
    ### (2) plot local pH for different setups ###
    dat = sample_local_pH_convection()
    plot_local_pH_convection('FigS10_local_pH_convection', dat)
    ##############################################


    #############################################################
    ### (3) plot surface concentrations of H2CCO, OH, CO with ###
    ###     varying pH and roughness                          ###
    csurf_kt, csurf_oh, csurf_co = get_surf_conc()
    plot_surf_concentrations('Fig4cde_surface_concentrations', csurf_co, csurf_oh, csurf_kt)
    #############################################################

    
