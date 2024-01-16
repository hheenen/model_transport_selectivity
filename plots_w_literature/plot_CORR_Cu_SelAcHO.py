#!/usr/bin/env python

from header import *

def pltargs(k):
    #clsk = {"1":clrs['lightblue'], "5":clrs['lightblue'], "87":clrs['orange'], "390":clrs['darkgray']}
    clsk = {"1":clrs['darkblue'], "5":clrs['lightblue'], "87":clrs['azurblue'], "390":clrs['lightblue']}
    if k.find('model') != -1:
        ls='-'; marker='None'; zorder=1
    else:
        ls='--'; marker='x'; zorder=99
    for t in clsk:
        if k.find(t) != -1:
            clr = clsk[t]
    return dict(ls=ls, marker=marker, color=clr, zorder=zorder)
   

if __name__ == "__main__":

    # TODO:
    # (x) reaction as header for plots
    # (x) refactor/document script: function to compute, etc.
    # (x) double check whehter percent in experimental output <-- mach mal!
    # (x) this plot change, color for roughness (blue-lightblue and fill for 1-5 oder colors for 87 and 390)
    

    #################################
    ### Parameters for simulation ###
    #################################

    ddG_des = 0.165
    dG_red = 1.1
    # dG --> G_des - G_ads 

    # Acetaldehyde on Cu
    cD_Acdh = 13.75e-10 # m2/s (10.3390/atmos11101057 average of two)
    cH_Acdh = 0.0759    # atm /mol (10.1016/S0045-6535(00)00505-1)
    dG_Acdh = [-0.6]    # 10.1002/anie.201508851
    dG_Acdh = -0.14 #-0.12, -0.15, -0.16 
    ddGred_Acdh = 0.96 #0.92
    ddG_des_Acdh = 0.165 # 0.1, 0.16, 0.165

    # effectively used model parameters
    cdict = {'Acdh_1':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_3':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_5':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_10':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_87':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             'Acdh_390':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh, 'ddG':ddG_des_Acdh, 'ddGred':ddGred_Acdh}, # Acdh for reasonable adsorption rates
             #'Acdh':  {'D':cD_Acdh, 'H':cH_Acdh, 'dG':dG_Acdh[0], 'ddG':ddG_des, 'ddGred':dG_red},
        }

    # calculation instructions
    sdict = {'Acdh_1':{"Us":np.arange(-1.35,-1.0,0.05), "rghs":[1], "mdls":[1]},
             'Acdh_3':{"Us":np.arange(-1.35,-1.0,0.05), "rghs":[3], "mdls":[1]},
             'Acdh_5':{"Us":np.arange(-1.35,-1.0,0.05), "rghs":[5], "mdls":[1]},
             #'Acdh_10':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[10], "mdls":[1]},
             #'Acdh_1':{"Us":np.arange(-1.3,-1.0,0.05), "rghs":[1], "mdls":[1]},
             'Acdh_87':{"Us":np.arange(-1.35,-1.0,0.05), "rghs":[87], "mdls":[1]},
             'Acdh_390':{"Us":np.arange(-1.35,-1.0,0.05), "rghs":[390], "mdls":[1]},
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
        
        # save potential vs. sel 
        
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        out_plot.update({k:np.array([dat[:,4], sel]).T})
    

    ##########################################
    ### Plot simulation against experiment ###
    ##########################################
    
    # prep experimental data -- from paper collection
    dat = load_Acdh_Cu_data()
    dat.update({'pc-Cu $\\rho$=1':dat.pop('COR@pc-Cu $\\rho$=1.0')})
    dat.update({'OD-Cu $\\rho$=87':dat.pop('COR@OD-Cu $\\rho$=87')})
    dat.update({'Cu-flower $\\rho$=390':dat.pop('COR@Cu-Flower $\\rho$=390')})
    dat.update({r'model $\rho=1$':out_plot['Acdh_1']})
    #dat.update({r'model $\rho=5$':out_plot['Acdh_5']})
    dat.update({r'model $\rho=87$':out_plot['Acdh_87']})
    dat.update({r'model $\rho=390$':out_plot['Acdh_390']})
    
    for k in dat:
        dat[k][:,1] *= 100.

    ls_args = {}
    for k in dat:
        ls_args.update({k:pltargs(k)})

    #ls_args['model $\\rho=1$']['color'] = 'w'
    #ls_args['model $\\rho=5$']['color'] = 'w'

   #ks = list(dat.keys())
   #cls = plt.cm.jet(np.linspace(0,1,len(ks)*2)[1::2])
   #ls_args = {ks[i]:dict(ls='--', marker='x', color=cls[i]) for i in range(len(ks))}
   #ls_args.update({r'model $\rho=1$':dict(ls='-', color=cls[0])})
   #ls_args.update({r'model $\rho=87$':dict(ls='-', color=cls[1])})
   #ls_args.update({r'model $\rho=390$':dict(ls='-', color=cls[2])})
   #ls_args = {ks[i]:dict(ls='--', marker='x', color='k') for i in range(len(ks))}
   #ls_args.update({r'model $\rho=1$':dict(ls='-', color='r')})
   #ls_args.update({r'model $\rho=5$':dict(ls='-', color='r')})
   #ls_args.update({r'model $\rho=87$':dict(ls='-', color='r')})
   #ls_args.update({r'model $\rho=390$':dict(ls='-', color='r')})
    xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity MeCHO / C$_{2\!+} \!(\%)$'

    extraleg = {k:ls_args[k] for k in ls_args if k.find('model') == -1}
    extraleg.update({'model':dict(ls='-', color='k')})
    
    filename = "ref_CORR_Cu_SelAcHO_pot"
    set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05, lrbt=[0.16,0.98,0.18,0.92],fsize=10)
    ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, \
        legpos=2, title='', extra_leg=extraleg)
        #legpos=2, title=r'CORR on Cu', extra_leg=extraleg)
    
    #ax.annotate("", xy=(2.5, 10), xytext=(2.5, 26), arrowprops={"arrowstyle":'<->', "color":'darkgray'})
   #ax.fill_between(x=dat[r'model $\rho=1$'][:,0], y1=dat[r'model $\rho=1$'][:,1], \
   #    y2=dat[r'model $\rho=5$'][:,1], fc=clrs['lightblue'], ec='w', alpha=0.3, zorder=0)
    #ax.annotate(r"$\rho=1-5$", xy=(-1.23, 27), size=7, color=clrs['lightblue'])
    #ax.set_ylim(ax.get_ylim()[0], 62)
    
    ax.annotate(r'$\rho$=1', xytext=(-1.31, 22), xy=(-1.26,20), color=clrs['darkblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['darkblue']))
    ax.annotate(r'$\rho$=87', xytext=(-1.23, 4), xy=(-1.15, 2), color=clrs['azurblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['azurblue']))
    ax.annotate(r'$\rho$=390', xytext=(-1.05, -1), xy=(-1.01,8), color=clrs['lightblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['lightblue']))
    #clsk = {"1":clrs['darkblue'], "5":clrs['lightblue'], "87":clrs['azurblue'], "390":clrs['lightblue']}
    
    add_sketch(ax, r'CO', r'C$_{2\!+}$', r'MeCHO', dxy=(0.12, 0.05))
    
    ax.annotate('(b)', xy=(0.08, 0.95), xycoords='figure fraction')
    ax.set_ylim(-4.87, 62) #hardcoded
    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=True) 


