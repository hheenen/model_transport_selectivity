#!/usr/bin/env python

"""

This script runs model for fig 2b in manuscript:
Methanol oxidation

"""

from header import *


def load_literature_MeOHox_Pt_data():
    """
      helper-function to load MeOHox data

    """
    dat_n = np.loadtxt("../literature_data/dat_H2CO_Pt_neg.txt")
    dat_p = np.loadtxt("../literature_data/dat_H2CO_Pt_pos.txt")

    #### H2CO on Pt oxidation
    selH2COn, selH2COp = dat_n[:,:2], dat_p[:,:2]
    selH2COn[:,0] /= selH2COn[:,0].min() # relative loading / roughness
    selH2COp[:,0] /= selH2COp[:,0].min() # relative loading / roughness

    return selH2COn, selH2COp
    

def plot_MeOHox_Pt(filename, dat, ls_args, sketch):
    """
      helper-function to plot MeOHox data

    """
    set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05,\
               lrbt=[0.16,0.98,0.185,0.92],fsize=10)
    
    xlabel = r'roughness $\rho$'; ylabel = r'selectivity CH$_2$O (%)'
    
    ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, \
        tag="", line=[], lsize=7)
    
    ax.set_xlim(0.7749999999999999, 5.7250000000000005) # hardcoded if model is not included
    #ax.set_ylim(21.13366790022691, 64.57458724284633) # hardcoded -- no model no neg scan
    
    ax.annotate('(b)', xy=(0.04, 0.95), xycoords='figure fraction')
    
    if sketch:
        #add_sketch(ax, r'MeOH', r'HCOOH + CO$_2$', r'H$_2$CO', dxy=(-0.3,-0.45))
        add_sketch(ax, r'MeOH', r'HCOOH+CO$_2$', r'CH$_2$O', dxy=(0.01,0.15))

    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=False) 


def make_plots_MeOHox_Pt():
    """
      helper-function to make full plots

    """

    # run model for MeOH ox on Pt example
    out_plot = run_model_for_example({'Fmdh':cdict}, {'Fmdh':sdict})
    out_plot['Fmdh'][:,0] /= out_plot['Fmdh'][:,0].min() # relative roughness
    out_plot['Fmdh'][:,1] *= 100.                        # convert percentage

    # load literature data
    selH2COn, selH2COp = load_literature_MeOHox_Pt_data()

    ls_args = {r'neg scan':dict(ls='--', marker='x', color='gray'), 
               r'pos scan':dict(ls='--', marker='x', color='k'), 
               r'Pt-NP':dict(ls='--', marker='x', color='k'), 
               r'model':dict(ls='-', color='r'),
               }

    # plot Fig2b
    filename = "Fig2b_MeOHox_Pt_SelFmdh"
    dat = {r'Pt-NP': selH2COp, 
        r'model':out_plot['Fmdh']}
    plot_MeOHox_Pt(filename, dat, ls_args, sketch=True)
    
    # plot FigS3b
    filename = "FigS3b_MeOHox_Pt_SelFmdh"
    dat = {r'neg scan': selH2COn,
        r'pos scan': selH2COp,
        r'model':out_plot['Fmdh']}
    plot_MeOHox_Pt(filename, dat, ls_args, sketch=False)


#################################
### Parameters for simulation ###
#################################

# dG --> G_des - G_ads 
ddG_des = 0.165      # guess
dG_red = 1.1         # guess

# Formaldehyde vs. HCOOH & CO2 --> digitize data
cD_Fmdh = 17.541e-10   # m2/s (10.1021/ie201944h0) fits to average of 10.3390/atmos11101057 
cH_Fmdh = 1.99e-4      # atm / M (10.1016/j.atmosenv.2010.05.044)
dG_Fmdh = [0.0]       # guess

# used model parameters
cdict = {'cD':[cD_Fmdh], 'cH':[cH_Fmdh], 'dG':[dG_Fmdh[0]], \
    'ddG_des':[ddG_des], 'ddG_red':[dG_red]}

# calculation parameters
sdict = {"Us":[-1.60], "rghs":np.arange(0.1,0.6,0.05), "mdls":[1]}


if __name__ == "__main__":
    make_plots_MeOHox_Pt()

