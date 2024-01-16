#!/usr/bin/env python

from header import *


if __name__ == "__main__":

    #################################
    ### Parameters for simulation ###
    #################################

    # dG --> G_des - G_ads 
    ddG_des = 0.165
    dG_red = 1.1

    # Formaldehyde vs. HCOOH & CO2 --> digitize data
    cD_Fmdh = 17.541e-10   # m2/s (10.1021/ie201944h0) fits to average of 10.3390/atmos11101057 
    cH_Fmdh = 1.99e-4      # atm / M (10.1016/j.atmosenv.2010.05.044)
    dG_Fmdh = [0.0]       # guess

    # effectively used model parameters
    cdict = {'Fmdh':  {'D':cD_Fmdh, 'H':cH_Fmdh, 'dG':dG_Fmdh[0], 'ddG':ddG_des, 'ddGred':dG_red},
        }

    # calculation instructions
    sdict = {'Fmdh':{"Us":[-1.60], "rghs":np.arange(0.1,0.6,0.05), "mdls":[1]},
        }
    

    ########################
    ### Simulation block ###
    ########################

    out_plot = {}
    for k in sdict:
        dGdes , dGads = get_adsdes_engs(cdict[k])

        print(k, dGdes, dGads)

        datafile = "model_data_examples_%s.pkl"%k
        # i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
        dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=[cdict[k]['ddGred']], \
            Ds=[cdict[k]['D']], Lxs=[50e-4], **sdict[k])
        
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        out_plot.update({k:np.array([dat[:,7], sel]).T})
    

    ##########################################
    ### Plot simulation against experiment ###
    ##########################################

    #### H2CO on Pt oxidation
    selH2COn, selH2COp = load_H2CO_Pt_data()
    selH2COn[:,0] /= selH2COn[:,0].min() # relative loading / roughness
    selH2COp[:,0] /= selH2COp[:,0].min() # relative loading / roughness
    out_plot['Fmdh'][:,0] /= out_plot['Fmdh'][:,0].min()
    out_plot['Fmdh'][:,1] *= 100.
        
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


