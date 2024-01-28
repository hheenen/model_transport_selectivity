#!/usr/bin/env python

from header import *
            


if __name__ == "__main__":

    # TODO:
    # (#) reaction as header for plots
    # (#) refactor script: function to compute, etc.
    # ( ) plot H2O2 needs the experimental data plotted in inset (to show cool effect)
    #     for this plot_xy needs to be giving out figure!


    #################################
    ### Parameters for simulation ###
    #################################

    # dG --> G_des - G_ads 
    ddG_des = 0.165
    dG_red = 1.1

    # H2O2 on Pt (Au check Sara's paper)
    cD_h2o2 = 15.6e-10     # m2/s --> 10.1039/C2CP40616K various values electrolyte depended
    cH_h2o2 = 1.20e-05     # atm / M  --> 10.1021/jp951168n
    dG_h2o2_Pt = [-0.25, -0.17]  # 10.1021/acs.jpcc.0c02127
    dG_h2o2_Au = [-0.34, 0.05]   # 10.1021/acs.jpcc.0c02127
    ddG_h2o2_des = 0.2 #0.35               # 10.1021/acs.jpcc.0c02127

    # effectively used model parameters
    cdict = {
             'H2O2_Pt':  {'cD':[cD_h2o2], 'cH':[cH_h2o2], 'dG':0.1, 'ddG_des':[ddG_h2o2_des], 'ddG_red':[dG_red]},
             #'H2O2_Au_lo':  {'cD':cD_h2o2, 'cH':cH_h2o2, 'dG':dG_h2o2_Au[1], 'ddG_des':ddG_h2o2_des, 'ddG_red':dG_red},
             #'H2O2_Au_ho':  {'cD':cD_h2o2, 'cH':cH_h2o2, 'dG':dG_h2o2_Au[0], 'ddG_des':ddG_h2o2_des, 'ddG_red':dG_red},
        }

    # calculation instructions
    sdict = {
             #'H2O2_Pt':{"Us":[-1.8], "rghs":[0.1,0.22,1.0], "mdls":[1]},
             'H2O2_Pt':{"Us":[-1.8], "rghs":np.arange(0.1,1.0,0.05), "mdls":[1]},
        }
    

    ########################
    ### Simulation block ###
    ########################

    out_plot = {}
    for k in sdict:
        dGdes , dGads = get_adsdes_engs(cdict[k])

        print(k, dGdes, dGads)

        datafile = "model_data_examples_%s.pkl"%k
        # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
        dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=cdict[k]['ddG_red'], \
            Ds=cdict[k]['cD'], Lxs=[50e-4], **sdict[k])
        # save roughness vs. sel 
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        out_plot.update({k:np.array([dat[:,7], sel]).T})
    

    ##########################################
    ### Plot simulation against experiment ###
    ##########################################

    #### ORR on Pt reduction
    selORR = load_ORR_Pt_data(); ks = list(selORR.keys())
    dat_all = {r'%.2f V vs. RHE'%k:selORR[k] for k in selORR}
    dat_all.update({r'pc-Pt': dat_all[r'0.58 V vs. RHE']})
   #if True:
   #    dat[r'Pt-disc'] = dat[r'0.58 V vs. RHE']
   #    del dat[r'0.53 V vs. RHE']
   #    del dat[r'0.48 V vs. RHE']
   #    del dat[r'0.43 V vs. RHE']
   #    del dat[r'0.58 V vs. RHE'] 
    dat_all.update({r'model':out_plot['H2O2_Pt']})
    dat_all[r'model'][:,1] *= 100.

    cls = plt.cm.jet(np.linspace(0,1,len(ks)))
    #ls_args = {r'ORR@Cu $^{[]}$':dict(ls='--', marker='x', color='k'), 
    ls_args = {r'%.2f V vs. RHE'%ks[i]:dict(ls='--', marker='x', color=cls[i], alpha=0.3) for i in range(len(ks))}
    ls_args[r'0.58 V vs. RHE']['alpha'] = 1.0
    ls_args.update({r'model':dict(ls='-', color='r')})
    ls_args.update({r'pc-Pt':dict(ls='--', marker='x', color='k')})
    
    dcoll = {'':[r'pc-Pt',r'model'], '_SI':[r'%.2f V vs. RHE'%k for k in selORR]+[r'model']}

    for tag in ['', '_SI']:
        dat = {k:dat_all[k] for k in dcoll[tag]}
        
        xlabel = r'roughness $\rho$'; ylabel = r'selectivity H$_2$O$_2$ (%)'
        
        set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05,\
                   lrbt=[0.16,0.98,0.185,0.92],fsize=10)
        
        ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7)#, title=r'ORR on Pt')
        ax.set_xlim(0.055, 1.045) # hardcoded if model is not included
        ax.set_ylim(3.7159276980511073, 35.53675756282097) # hardcoded -- no model
        ax.annotate('(a)', xy=(0.04, 0.95), xycoords='figure fraction')

        if tag == '':
            add_sketch(ax, r'O$_2$', r'H$_2$O', r'H$_2$O$_2$')
 
        filename = "ref_ORR_Pt_SelH2O2"+tag
        writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=True) 
    
