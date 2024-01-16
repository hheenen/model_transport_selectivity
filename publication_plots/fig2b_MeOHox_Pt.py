#!/usr/bin/env python

"""

This script runs model for fig 2b in manuscript:
Methanol oxidation

"""

from header import *

    
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


#   def make_plot_MeOHox_Pt():
#   """
#     helper-function to make full plot

#   """
#   pass

    

    ########################
    ### Simulation block ###
    ########################

    out_plot = run_model_for_example({'Fmdh':cdict}, {'Fmdh':sdict})
    out_plot['Fmdh'][:,0] /= out_plot['Fmdh'][:,0].min() # relative roughness
    out_plot['Fmdh'][:,1] *= 100.                        # convert percentage

    ##########################################
    ### Plot simulation against experiment ###
    ##########################################

    #### H2CO on Pt oxidation
    selH2COn, selH2COp = load_H2CO_Pt_data()
    selH2COn[:,0] /= selH2COn[:,0].min() # relative loading / roughness
    selH2COp[:,0] /= selH2COp[:,0].min() # relative loading / roughness

    assert False
    # ( ) sort out central data here!
    # ( ) make plot function and have central

    ls_args = {r'neg scan':dict(ls='--', marker='x', color='gray'), 
               r'pos scan':dict(ls='--', marker='x', color='k'), 
               r'Pt-NP':dict(ls='--', marker='x', color='k'), 
               r'model':dict(ls='-', color='r'),
               }
    dat_all = {r'neg scan':selH2COn,
           r'pos scan':selH2COp,
           r'Pt-NP':selH2COp,
           r'model':out_plot['Fmdh']
           }
    dcoll = {'':[r'Pt-NP',r'model'], '_SI':[r'neg scan', r'pos scan', r'model']}

    for tag in ['', '_SI']:
        dat = {k:dat_all[k] for k in dcoll[tag]}

        xlabel = r'roughness $\rho$'; ylabel = r'selectivity CH$_2$O (%)'
        
        set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05,\
                   lrbt=[0.16,0.98,0.185,0.92],fsize=10)
        
        ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7)#, title=r'MeOH oxidation on Pt')
        ax.set_xlim(0.7749999999999999, 5.7250000000000005) # hardcoded if model is not included
        #ax.set_ylim(21.13366790022691, 64.57458724284633) # hardcoded -- no model no neg scan
        ax.annotate('(b)', xy=(0.04, 0.95), xycoords='figure fraction')
        
        if tag == '':
            #add_sketch(ax, r'MeOH', r'HCOOH + CO$_2$', r'H$_2$CO', dxy=(-0.3,-0.45))
            add_sketch(ax, r'MeOH', r'HCOOH+CO$_2$', r'CH$_2$O', dxy=(0.01,0.15))
 
        filename = "ref_MeOHox_Pt_SelFmdh"+tag
        writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=True) 


