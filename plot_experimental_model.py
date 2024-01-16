#!/usr/bin/env python

"""

plotting some main data

"""

import os, sys
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from sim_diffusion_selectivity import sample_data

sys.path.append("experimental_reference")
from read_data import *
from model_plots import plot_xy


def _get_adsdes_engs(cd):
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


if __name__ == "__main__":

    
    ################################
    #### Literature and adjusted model parameters

    ddG_des = 0.165
    dG_red = 1.1
    # dG --> G_des - G_ads 

    # CO on Cu
    cH_CO = 1100           # Henry's constant for CO atm/M  (10.5194/acp-15-4399-2015)
    cD_CO = 20.3e-10       # m2/s              (Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems)
    dG_CO_111 = [0.1, 0.28] # (10.1039/C0EE00071J , 10.1021/acsenergylett.0c01751)
    dG_CO_211 = [0.4]       # 10.1021/acsenergylett.0c01751
    dG_red_CO = 1.45

    # Acetaldehyde on Cu
    cD_Acdh = 13.75e-10 # m2/s (10.3390/atmos11101057 average of two)
    cH_Acdh = 0.0759    # atm /mol (10.1016/S0045-6535(00)00505-1)
    dG_Acdh = [-0.6]    # 10.1002/anie.201508851


    # H2O2 on Pt (Au check Sara's paper)
    cD_h2o2 = 15.6e-10     # m2/s --> 10.1039/C2CP40616K various values electrolyte depended
    cH_h2o2 = 1.20e-05     # atm / M  --> 10.1021/jp951168n
    dG_h2o2_Pt = [-0.25, -0.17]  # 10.1021/acs.jpcc.0c02127
    dG_h2o2_Au = [-0.34, 0.05]   # 10.1021/acs.jpcc.0c02127
    ddG_h2o2_des = 0.2 #0.35               # 10.1021/acs.jpcc.0c02127
    

    # Formaldehyde vs. HCOOH & CO2 --> digitize data
    cD_Fmdh = 17.541e-10   # m2/s (10.1021/ie201944h0) fits to average of 10.3390/atmos11101057 
    cH_Fmdh = 1.99e-4      # atm / M (10.1016/j.atmosenv.2010.05.044)
    dG_Fmdh = [0.0]       # guess

    dG_Acdh = -0.14 #-0.12, -0.15, -0.16 
    ddGred_Acdh = 0.96 #0.92
    ddG_des_Acdh = 0.165 # 0.1, 0.16, 0.165

    ### effectively used model parameters
    cdict = {'CO_111':{'D':cD_CO, 'H':cH_CO, 'dG':0.1, 'ddG':ddG_des, 'ddGred':dG_red_CO},
             #'CO_211':{'D':cD_CO, 'H':cH_CO, 'dG':0.4, 'ddG':ddG_des, 'ddGred':dG_red_CO},
             'Acdh_1':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_3':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_5':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_10':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_87':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_390':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             #'Acdh':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh[0], 'ddG':ddG_des, 'ddGred':dG_red},
             'H2O2_Pt':  {'D':cD_h2o2, 'H':cH_h2o2, 'dG':0.1, 'ddG':ddG_h2o2_des, 'ddGred':dG_red},
             #'H2O2_Au_lo':  {'D':cD_h2o2, 'H':cH_h2o2, 'dG':dG_h2o2_Au[1], 'ddG':ddG_h2o2_des, 'ddGred':dG_red},
             #'H2O2_Au_ho':  {'D':cD_h2o2, 'H':cH_h2o2, 'dG':dG_h2o2_Au[0], 'ddG':ddG_h2o2_des, 'ddGred':dG_red},
             'Fmdh':  {'D':cD_Fmdh, 'H':cH_Fmdh, 'dG':dG_Fmdh[0], 'ddG':ddG_des, 'ddGred':dG_red},
        }

    sdict = {'CO_111':{"Us":[-1.26], "rghs":np.arange(1,5,0.2), "mdls":[1]}, # test CO
    #sdict = {'CO_111':{"Us":[-1.26], "rghs":[1,3,5], "mdls":[1]}, # test CO
             'Acdh_1':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[1], "mdls":[1]},
             'Acdh_3':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[3], "mdls":[1]},
             'Acdh_5':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[5], "mdls":[1]},
             #'Acdh_10':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[10], "mdls":[1]},
             #'Acdh_1':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[1], "mdls":[1]},
             'Acdh_87':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[87], "mdls":[1]},
             'Acdh_390':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[390], "mdls":[1]},
             #'H2O2_Pt':{"Us":[-1.8], "rghs":[0.1,0.22,1.0], "mdls":[1]},
             'H2O2_Pt':{"Us":[-1.8], "rghs":np.arange(0.1,1.0,0.05), "mdls":[1]},
             'Fmdh':{"Us":[-1.60], "rghs":np.arange(0.1,0.6,0.05), "mdls":[1]},
        }

    out_plot = {}
    for k in sdict:
    #for k in ['Acdh']:
        dGdes , dGads = _get_adsdes_engs(cdict[k])

        print(k, dGdes, dGads)

        datafile = "model_data_examples_%s.pkl"%k
        dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=[cdict[k]['ddGred']], \
            Ds=[cdict[k]['D']], Lxs=[50e-4], **sdict[k])
        # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
# 
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        #if k == 'Acdh':
        if k.find('Acdh_') != -1:
            out_plot.update({k:np.array([dat[:,4], sel]).T})
        elif k != 'Acdh':
            out_plot.update({k:np.array([dat[:,7], sel]).T})
        print(sel)
#  #    print(sel[:3] / sel[3:6])
#  #    print(sel[:3] / sel[6:9])

    ################################
    #### Experimental example data
    

    #### CO reduction
##  selCOr, selCOir = load_CO_Cu_data()
##  av = np.mean([selCOr[:,0], selCOir[:,0]], axis=0)
##  err = np.zeros(selCOr.shape); err[:,0] = np.std([selCOr[:,0], selCOir[:,0]], axis=0)
##  selCO = np.array([av, selCOr[:,1]]).T

##  #dat = {r'CO$_2$RR@Cu $^{[]}$':selCO,
##  dat = {r'CO$_2$RR@Cu':selCO,
##         #r'CO@Cu $^{[]}$ - alt':selCOir,
##         r'model':out_plot['CO_111'],
##         }
##  err = {r'CO$_2$RR@Cu $^{[]}$':err}
##  for k in dat:
##      dat[k][:,1] *= 100.

##  ls_args = {r'CO$_2$RR@Cu':dict(ls='--', marker='x', color='k'), 
##  #ls_args = {r'CO$_2$RR@Cu $^{[]}$':dict(ls='--', marker='x', color='k'), 
##             #r'CO$_2$RR@Cu $^{[]}$ - alt':dict(ls=':', marker='x', color='k'),
##             r'model':dict(ls='-', color='r'),
##             }
##  xlabel = 'roughness'; ylabel = r'selectivity CO (%)'
##  
##  plot_xy("ref_roughness_selectivity_CO", dat, ylabel, xlabel, ls_args, err, tag="", line=[], fsize=10, lsize=7)
    # NOTE: add arrow to say that sites may play a role as well

    
    #### ORR on Pt reduction
##  selORR = load_ORR_Pt_data(); ks = list(selORR.keys())
##  dat = {r'%.2f V vs. RHE'%k:selORR[k] for k in selORR}
##  dat[r'ORR@Pt'] = dat[r'0.58 V vs. RHE']
##  del dat[r'0.53 V vs. RHE']
##  del dat[r'0.48 V vs. RHE']
##  del dat[r'0.43 V vs. RHE']
##  del dat[r'0.58 V vs. RHE'] 
##  dat.update({r'model':out_plot['H2O2_Pt']})
##  dat[r'model'][:,1] *= 100.

##  cls = plt.cm.jet(np.linspace(0,1,len(ks)))
##  #ls_args = {r'ORR@Cu $^{[]}$':dict(ls='--', marker='x', color='k'), 
##  ls_args = {r'%.2f V vs. RHE'%ks[i]:dict(ls='--', marker='x', color=cls[i]) for i in range(len(ks))}
##  ls_args.update({r'model':dict(ls='-', color='r')})
##  ls_args.update({r'ORR@Pt':dict(ls='--', marker='x', color='k')})

##  xlabel = 'roughness'; ylabel = r'selectivity H$_2$O$_2$ (%)'
##  
##  plot_xy("ref_roughness_selectivity_H2O2", dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7)
    #plot_xy("ref_roughness_selectivity_H2O2", dat, ylabel, xlabel, ls_args, tag="ORR@Pt $^{[]}$", line=[], fsize=10, lsize=7)
    
    # ( ) same plots for H2O2 and Fmdh (for H2O2 more potentials, for Fmdh relative roughness)
    
    
##  #### H2CO on Pt oxidation
##  selH2COn, selH2COp = load_H2CO_Pt_data()
##  selH2COn[:,0] /= selH2COn[:,0].min() # relative loading / roughness
##  selH2COp[:,0] /= selH2COp[:,0].min() # relative loading / roughness
##  out_plot['Fmdh'][:,0] /= out_plot['Fmdh'][:,0].min()
##  dat = {#r'MeOHox@Pt neg $^{[]}$':selH2COn,
##         #r'MeOHox@Pt pos $^{[]}$':selH2COp,
##         r'MeOHox@Pt':selH2COp,
##         r'model':out_plot['Fmdh']
##         }
##  dat[r'model'][:,1] *= 100.
##  ls_args = {r'MeOHox@Pt neg $^{[]}$':dict(ls='--', marker='x', color='gray'), 
##             r'MeOHox@Pt pos $^{[]}$':dict(ls='--', marker='x', color='k'), 
##             r'MeOHox@Pt':dict(ls='--', marker='x', color='k'), 
##             r'model':dict(ls='-', color='r'),
##             }
##  xlabel = 'mass-normalized roughness'; ylabel = r'selectivity H$_2$CO (%)'
##  
##  plot_xy("ref_roughness_selectivity_Fmdh", dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7)
    

    #### AcOH on Cu
    # TODO: double check whehter percent in experimental output <-- mach mal!
    # TODO: this plot change, color for roughness (blue-lightblue and fill for 1-5 oder colors for 87 and 390)
    # TODO: plot H2O2 needs the experimental data plotted in inset (to show cool effect)
    # TODO: for this plot_xy needs to be giving out figure!
    # TODO: split this script into a script PER plot
    # TODO: Look at MeOH ox paper again (pos/neg scan)
    # TODO: reaction as header for plots
##  dat = load_Acdh_Cu_data()
##  dat.update({r'model $\rho=1$':out_plot['Acdh_1']})
##  dat.update({r'model $\rho=5$':out_plot['Acdh_5']})
##  dat.update({r'model $\rho=87$':out_plot['Acdh_87']})
##  dat.update({r'model $\rho=390$':out_plot['Acdh_390']})
##  
##  for k in dat:
##      dat[k][:,1] *= 100.
##  #dat[r'model'][:,1] *= 100.
##  
##  ks = list(dat.keys())
##  cls = plt.cm.jet(np.linspace(0,1,len(ks)*2)[1::2])
## #ls_args = {ks[i]:dict(ls='--', marker='x', color=cls[i]) for i in range(len(ks))}
## #ls_args.update({r'model $\rho=1$':dict(ls='-', color=cls[0])})
## #ls_args.update({r'model $\rho=87$':dict(ls='-', color=cls[1])})
## #ls_args.update({r'model $\rho=390$':dict(ls='-', color=cls[2])})
##  ls_args = {ks[i]:dict(ls='--', marker='x', color='k') for i in range(len(ks))}
##  ls_args.update({r'model $\rho=1$':dict(ls='-', color='r')})
##  ls_args.update({r'model $\rho=5$':dict(ls='-', color='r')})
##  ls_args.update({r'model $\rho=87$':dict(ls='-', color='r')})
##  ls_args.update({r'model $\rho=390$':dict(ls='-', color='r')})
##  xlabel = r'$\rm U_{SHE}~(V)$'; ylabel = r'selectivity MeCHO (%)'

##  extraleg = {'COR@Cu':dict(ls='--', color='k'), 'model':dict(ls='-', color='r')}
##  
##  plot_xy("ref_pot_selectivity_Acdh", dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, legpos=2, extra_leg=extraleg)


    #### Ac on CuPd
    # load acetate model
    
    # load experimental reference # TODO: save rho_A = x and rho_I = y plot that as the cool U-plot's then you think about model!
    dat = load_Ac_CuPd_data()
    print(dat)


