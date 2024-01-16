#!/usr/bin/env python

"""

This script only produces 1D cuts for each parameter

"""

import os, sys
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
sys.path.append("../.")
from sim_diffusion_selectivity import sample_data

from model_plots import plot_xy, plot_xy_ax, set_plotting_env, writefig, golden_ratio, clrs

def get_adsdes_engs(cd):
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


def adjust_param_for_plot(p, param):
    if p == 'dG':
        param = np.array(param) * -1 # let's call it adsorption free energy (previously was desorption
    elif p == 'ddG_red':
        param = np.array(param) * param_std['U']*0.5
    elif p == 'Lx':
        param = np.array(param) * 1e4
    elif p == 'cD':
        param = np.array(param) * 1e10
    return param


def plot_model_parameter_selectivity_1D(filename, dat, pstd, datm2={}):
    set_plotting_env(width=3.37*2,height=3.37/ golden_ratio *1.05,\
               lrbt=[0.08,0.98,0.19,0.98],fsize=9.0)
        
    # output dat: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1

    # rgh, Dx, Lx, Gads, dGdes, dGred*U*0.5
    klabel = {'dG':r'$\Delta G _{\mathrm{ads}}$ (eV)', 
              'ddG_des':r'$\Delta G ^{\ddag} _{\mathrm{des}}$ (eV)',
              'ddG_red':r'$\Delta G ^{\ddag} _{\mathrm{red}} (U, \alpha)$ (eV)', 
              'cD':r'$D_{\mathrm{diff}}$ ($10^{10}$ $\frac{\mathrm{m}^2}{s}$)', #m$^2$ s$^{-1}$)', 
              'Lx':r'$L_{\mathrm{diff}}$ ($\mathrm{\mu}$m)', 
              'rgh':r'roughness $\rho$' }

    tag = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

    fig, ax = plt.subplots(1,6)

    for i, k in enumerate(['dG', 'ddG_des', 'ddG_red', 'cD', 'Lx', 'rgh']):
        d = dat[k]
        ax[i].plot(d[:,0], d[:,1], color='k')
        ax[i].set_xlabel(klabel[k])
        ax[i].set_ylim(-0.02*100., 1.02*100.)
        ax[i].plot(pstd[k], 0.514*100., marker='o', color='k', markersize=2)
        ax[i].annotate(tag[i], xy=(0.1,0.85), xycoords='axes fraction')
        if k in datm2: # other model
            dd = datm2[k]
            ax[i].plot(dd[:,0], dd[:,1], color='k', ls='--')
    
    [ax[i].set_yticklabels([]) for i in range(1,6)]
    plt.subplots_adjust(wspace=0.05)

    if len(datm2) > 0:
        ax[1].annotate(r'model b', xy=(0.7,0.8), xytext=(0.1,0.7), size=7, \
            xycoords='axes fraction', arrowprops={"arrowstyle":'-', "color":'k', 'lw':1.0})

    ax[0].set_ylabel(r'selectivity $C_g$ (%)')
    
    writefig(filename ,folder='output', write_pdf=False, write_png=True) 


    


if __name__ == "__main__":

    param_std = dict(
        cH = [1],  # Henry's constant = 1
        cD = [20.0e-10],    # m2/s
        Lx = [50e-4],       # cm (50e-6m)
        dG = [-0.25],
        ddG_des = [0.15],
        ddG_red = [0.815],
        rgh = [1],
        U = [-1.0],
        )
    
    # ranges - rough later refine in case there are sensitive regions with maximum effect
    ranges = dict(
        cD = np.arange(10,61,5)*1e-10, # typical D(aq), \cite{Cussler}
        Lx = np.arange(10,100,5)*1e-4,    # ?! need range from levich equation estimate
        dG = np.arange(-0.45, -0.05, 0.05),
        ddG_des = np.arange(0.1,0.5,0.05), # around 0.15 \cite{Heenen2022}
        ddG_red = np.arange(0.65,1.05,0.05),
        rgh = np.array([1,2,3,5,8,10,15,20,25,30,35,40,50,60]),
        )

    
    ########################
    ### Simulation block ###
    ########################

    out = {}
    for mdl in [1,2]: # include model dependence in parameter screening
        out_plot = {}
        for p in ['cD', 'Lx', 'dG', 'ddG_des', 'ddG_red', 'rgh']:
            if p not in ['dG', 'ddG_des']: # complex workaround
                dGdes , dGads = get_adsdes_engs(param_std)
                dGdes = [dGdes]; dGads = [dGads]
            else:
                dGdes = []; dGads = []
                for e in ranges[p]:
                    pstd = deepcopy(param_std); pstd.update({p:[e]})
                    ddGdes , ddGads = get_adsdes_engs(pstd)
                    dGdes.append(ddGdes); dGads.append(ddGads)
                    
 
            params = deepcopy(param_std)
            params.update({p:ranges[p]})
 
            datafile = "model_data_screen_%s.pkl"%p
            # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
            if p not in ['dG', 'ddG_des']:
                dat = sample_data(datafile, rdes=dGdes, rads=dGads, 
                    rred=params['ddG_red'], Ds=params['cD'], Lxs=params['Lx'], 
                    Us=params['U'], rghs=params['rgh'], mdls=[mdl])
            else: # complex workaround for matched values
                dat = np.zeros((len(params[p]), 16))
                for i in range(len(params[p])):
                    dat[i,:] = sample_data(datafile, rdes=[dGdes[i]], rads=[dGads[i]], 
                        rred=params['ddG_red'], Ds=params['cD'], Lxs=params['Lx'], 
                        Us=params['U'], rghs=params['rgh'], mdls=[mdl])[0].tolist()
 
            # save roughness vs. sel (in %)
            sel = (dat[:,13]/ dat[:,[13,14]].sum(axis=1)) * 100.
            print(p, params[p], sel)
 
 
            params[p] = adjust_param_for_plot(p, params[p])
            
            out_plot.update({p:np.array([params[p], sel]).T})
        
        out.update({mdl:out_plot})


    pstd_plot = {p:adjust_param_for_plot(p, deepcopy(param_std[p])) for p in param_std}

    # out[2] == direct reduction to C* ; out[1] == reduction and desorptoin to C_aq
    plot_model_parameter_selectivity_1D('model_parameter_sel_1D', out[2], pstd_plot, out[1])
  
    # specific parameter output for selectivity change from 13.75e-10 to 20.3e-10
 #  params = deepcopy(param_std)
 #  dat = sample_data(datafile, rdes=dGdes, rads=dGads, 
 #      rred=params['ddG_red'], Ds=[13.75e-10, 20.3e-10], Lxs=params['Lx'], 
 #      Us=params['U'], rghs=params['rgh'], mdls=[mdl])
 #  sel = (dat[:,13]/ dat[:,[13,14]].sum(axis=1)) * 100.
 #  print(sel)


