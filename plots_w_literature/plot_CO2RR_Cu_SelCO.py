#!/usr/bin/env python

from header import *


if __name__ == "__main__":

    # TODO:
    # (x) reaction as header for plots
    # ( ) line where site-energy reduces as well
    # (x) add arrow to say that sites may play a role as well
    # (x) refactor/document script: function to compute, etc.
 
    
    #################################
    ### Parameters for simulation ###
    #################################

    # dG --> G_des - G_ads (i.e. positive means adsorption is exothermic)
    ddG_des = 0.165

    # parameters for CO on Cu
    cH_CO = 1100           # Henry's constant for CO atm/M  (10.5194/acp-15-4399-2015)
    cD_CO = 20.3e-10       # m2/s              (Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems)
    dG_CO_111 = [0.1, 0.28] # (10.1039/C0EE00071J , 10.1021/acsenergylett.0c01751)
    dG_CO_211 = [0.4]       # 10.1021/acsenergylett.0c01751
    dG_red_CO = 1.45        # NOTE: changed for CO selectivity to 1.5

    # effectively used model parameters
    cdict = {
        'CO_111':{'cD':[cD_CO], 'cH':[cH_CO], 'dG':[0.1], 'ddG_des':[ddG_des], 'ddG_red':[dG_red_CO]},
        #'CO_211':{'cD':cD_CO, 'cH':cH_CO, 'dG':0.4, 'ddG_des':ddG_des, 'ddG_red':dG_red_CO},
        }

    # calculation instructions
    sdict = {
        'CO_111':{"Us":[-1.26], "rghs":np.arange(1,5,0.2), "mdls":[1]}, # test CO
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
    
 #  # prep experimental data -- evaluate error from capacitive and current based roughness
 #  selCOr, selCOir = load_CO_Cu_data()
 #  av = np.mean([selCOr[:,0], selCOir[:,0]], axis=0)
 #  err = np.zeros(selCOr.shape); err[:,0] = np.std([selCOr[:,0], selCOir[:,0]], axis=0)
 #  selCO = np.array([av, selCOr[:,1]]).T

 #  label_exp =  r'experiment$^{[]}$' #r'CO$_2$RR@Cu' # r'CO$_2$RR@Cu $^{[]}$'
 #  dat = {label_exp:selCO, r'model':out_plot['CO_111']}
 #  err = {label_exp:err}
 #  for k in dat:
 #      dat[k][:,1] *= 100.

 #  ls_args = {label_exp:dict(ls='--', marker='x', color='k'), 
 #             r'model':dict(ls='-', color='r'),
 #             }
 #  xlabel = 'roughness'; ylabel = r'selectivity CO (%)'
 #  
 #  filename = "ref_CO2RR_Cu_SelCO"
 #  set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.16,0.98,0.17,0.92],fsize=10)
 #  ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, err, tag="", line=[], fsize=10, lsize=7, title=r'CO$_2$RR on Cu')
 #  ax.annotate("", xy=(2.5, 7), xytext=(2.5, 26), arrowprops={"arrowstyle":'<->', "color":'darkgray'})
 #  ax.annotate(r"$\Delta =$ more steps?", xy=(2.6, 15), size=7, color='darkgray')
 #  writefig(filename ,folder='output', write_pdf=False, write_png=True) 

    
    ##########################################
    ### other plot potential dependence CO ###
    ##########################################

    ### prep experimental data
    selCO, rghCO, phCO = load_CO_Cu_data_collection()
    kout = ['Huang-110', 'Huang-100', 'Wang-NSph-Cu+pH14.0']
    selCO = {k:selCO[k] for k in selCO if k not in kout}
    rghCO['Huang-OD-Cu'] = 32 # Huang-OD-Cu is roughly roughness = 30 estimated from current; Huang_100 is also much higher in roughness...!
    rghCO['Kuhl-pc-Cu'] = 1.1
    rghCO['Huang-111'] = 1.2
    rghCO['Huang-110'] = 1.3 # roughness is actually 4 according to estimate
    rghCO['Huang-100'] = 1.4 # roughness is actually 4 according to estimate
    rghCO['Kanan-pc-Cu'] = 30.1 # from total current totally visible in paper that roughness = 30
    rghCO = {k:float(rghCO[k]) for k in rghCO}

    print(rghCO)

    clsk = {1:clrs['darkblue'], 1.1:clrs['azurblue'], 1.2:clrs['lightblue2'], 1.3:clrs['lightblue3'], 1.4:'b',
        30:'r', 30.1:clrs['orange'], 32:clrs['darkyellow'],
        198:clrs['darkred'], 475.0:clrs['darkred']} #"390":clrs['darkgray']}

    selCO = {k:selCO[k] for k in selCO if phCO[k] < 10.}
    #selCO = {k:selCO[k] for k in selCO if phCO[k] > 10.}
    rghCO = {k:rghCO[k] for k in selCO}
    phCO = {k:phCO[k] for k in selCO}
    ls_args = {k:dict(ls='--', marker='x', color=clsk[rghCO[k]]) for k in selCO}

  # del selCO['Huang-OD-Cu'] # for Vanessa build-up
  # del selCO['Kanan-pc-Cu'] # for Vanessa build-up


    ### simulated data
    # Important: adjusted alhpa to 0.4 in function "adjust_SHE_engs" file 
    # "model_COR_Ac/mkm/mkm_energetics.py", by hand for simplicity ! NOTE: not necessary anymore!
    for rgh in [1,30]:#,198]: #475]:
        #datafile = "model_data_examples_CO_pot_alpha_hand.pkl"
        datafile = "model_data_examples_CO_pot.pkl"
        # output: i_model, eng_des, eng_ads, eng_red, U_SHE, Dx, Lx, roughness, *, A, B, C, D, p1, p2, conc1
        dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=[1.5], #][cdict[k]['ddG_red']], \
        #dat = sample_data(datafile, rdes=[dGdes], rads=[dGads], rred=[1.37], \
            Ds=cdict[k]['cD'], Lxs=[50e-4], Us=np.arange(-1.6,-0.6,0.05), rghs=[rgh], mdls=[1])
        # save roughness vs. sel 
        sel = dat[:,13]/ dat[:,[13,14]].sum(axis=1)
        ind = (np.isnan(sel) == False)
        selCO.update({'sim_%i'%rgh:np.array([dat[ind,4], sel[ind]]).T})
        ls_args.update({'sim_%i'%rgh:dict(ls='-', lw=2.5, color=clsk[rgh])})

    extraleg = {r'Cu(111) $\rho=$1':dict(ls='--', marker='x', color=clsk[1.2]),
        #r'Cu(110) $\rho \approx$1':dict(ls='--', marker='x', color=clsk[1.3]),
        r'pc-Cu $\rho=$1':dict(ls='--', marker='x', color=clsk[1.1]),
       #r'model $\rho=$1':dict(ls='-', color=clsk[1]),
        #r'OD-Cu $\rho \approx$30 ':dict(ls='--', marker='x', color=clsk[30.1]),
        r'OD-Cu $\rho=$30':dict(ls='--', marker='x', color=clsk[30.1]),
        r'OD-Cu $\rho=$32':dict(ls='--', marker='x', color=clsk[32]),
     #  r'OD-Cu $\rho \approx$30':dict(ls='--', marker='x', color=clsk[30.1]),
     #  r'OD-Cu $\rho \approx$32':dict(ls='--', marker='x', color=clsk[32]),
       #r'model $\rho = $30':dict(ls='-', color=clsk[30]),
        r'model':dict(ls='-', color='r'),
    }
 ###ks = list(selCO.keys())
 ###cls = plt.cm.jet(np.linspace(0,1,len(ks)*2)[1::2])
 ###ls_args = {ks[i]:dict(ls='--', marker='x', color=cls[i]) for i in range(len(ks))}
 ###extraleg = []
 
    for k in selCO:
        selCO[k][:,1] *= 100.


    ### plotting
    filename = "ref_CO2RR_Cu_SelCO_pot"
    xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity CO / C$_{1} \! + \! \mathrm{C}_{2\!+} \! (\%)$'# + 2\!\!+}$ (%)'
    #xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity CO / C$_{1\!+\!2+}$ (%)'
    
    set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05, lrbt=[0.17,0.98,0.18,0.92],fsize=10)
    ax = plot_xy_ax(selCO, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, \
        legpos=4, title='', extra_leg=extraleg)
        #legpos=4, title=r'CO$_2$RR on Cu', extra_leg=extraleg)

    #ax.annotate('(a)', xy=(0.04, 0.95), xycoords='figure fraction')
    ax.annotate('(a)', xy=(0.08, 0.95), xycoords='figure fraction')
    ax.annotate(r'$\rho$=1', xytext=(-1.62, 25), xy=(-1.45,25), color=clsk[1], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clsk[1]))
    ax.annotate(r'$\rho$=30', xytext=(-1.25, 5), xy=(-1.35,10), color=clsk[30], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clsk[30]))
    
    #add_sketch(ax, r'CO$_2$', r'C$_{1\!+\!2+}$', r'  CO', dxy=(0.22,-0.2))
    #add_sketch2(ax, r'CO$_2$', r'C$_{1\!+\!2\!\!+}$', r'  CO') #dxy=(0.22,-0.2))
    #add_sketch2(ax, r'CO$_2$', r'C$_{1 + 2\!\!+}$', r'  CO') #dxy=(0.22,-0.2))
    add_sketch2(ax, r'CO$_2$', r'C$_{1} \! + \! \mathrm{C}_{2\!+}$', r'  CO') #dxy=(0.22,-0.2))

    # ax.set_xlim(ax.get_xlim()[0], -0.65)
    ax.set_xlim(-1.70252, -0.65) #hardcoded
    #ax.set_xlim(-1.78, -0.65) #hardcoded
    ax.set_ylim(-4.93, 104.99) #hardcoded
    #ax.yaxis.set_label_coords(-0.135, 0.42)
    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=True)

