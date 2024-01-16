#!/usr/bin/env python

from header import *

'Koper-Cu_flat+pH6.8', 
'Koper-OD-Cu_sphere+pH6.8', 
'Koper-OD-Cu_cubes+pH6.8', 
'Koper-OD-Cu_dendrites+pH6.8', 
'Koper-Cu90Zn10_flat+pH6.8', 
'Koper-OD-Cu90Zn10_sphere+pH6.8', 
'Koper-OD-Cu90Zn10_cubes+pH6.8', 
'Koper-OD-Cu90Zn10_cauliflower+pH6.8', 
'Koper-Cu75Zn25_flat+pH6.8', 
'Koper-OD-Cu75Zn25_sphere+pH6.8', 
'Koper-OD-Cu75Zn25_cubes+pH6.8', 
'Koper-OD-Cu75Zn22_cauliflower+pH6.8'



def kf(k):
    kk = k[6:].split('+')[0]
    return kk
    

if __name__ == "__main__":
    

    d_sel, d_rgh = get_CuZn_alloy_CO()
    ks = list(set([k[6:].split('_')[0].split('-')[-1] for k in d_sel]))
    cls = plt.cm.cool(np.linspace(0,1,len(ks)))
    clrs = {ks[i]:cls[i] for i in range(len(ks))}
    ylabel = r'selectivity CO / C$_{1+2}$ (%)'

    # (1) sel vs. rgh - sort by alloy
    iU = 2
    dat = {k:[[d_rgh[kk], d_sel[kk][iU,1]] for kk in d_sel if kk.split('-')[-1].split('_')[0] == k] for k in ks}
    for k in dat:
        dat[k] = np.array(dat[k])
        dat[k][:,1] *= 100.
    

    xlabel = r'roughness $\rho$'
    ls_args = {k:dict(ls='--', marker='x', markersize=7, markeredgewidth=2, color=clrs[k]) for k in ks}

    filename = "ref_CO2RR_CuZn_SelCO"
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.16,0.98,0.17,0.92],fsize=10)
    ax = plot_xy_ax(deepcopy(dat), ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=6, \
        legpos=1, title=r'')#, extra_leg=extraleg)
    ax.annotate("%.2f U$_{SHE}$"%d_sel['Koper-OD-Cu75Zn25_cauliflower+pH6.8'][iU,0], size=8, xy=(0.8, 0.95), xycoords='figure fraction')
    writefig(filename ,folder='output', write_pdf=False, write_png=True)


    # (2) sel vs. U
    dat = {r"%s $\rho=%.1f$"%(kf(k), d_rgh[k]):d_sel[k] for k in d_rgh}
    for k in dat:
        dat[k] = np.array(dat[k])
        dat[k][:,1] *= 100
    
    xlabel = r'$\rm U_{SHE}~(V)$';
    print(dat.keys())
    lss = {'flat':'-', 'sphere':'-.', 'cubes':'--', 'cauliflower':':', 'dendrites':':'}
    ls_args = {k:dict(ls=lss[k.split()[0].split('_')[-1]], marker='x', markersize=7, markeredgewidth=2, color=clrs[k.split('-')[-1].split('_')[0]]) for k in dat}
    
    filename = "ref_CO2RR_CuZn_SelCO_pot"
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.16,0.98,0.17,0.92],fsize=10)
    ax = plot_xy_ax(deepcopy(dat), ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=6, \
        legpos=1, title=r'')#, extra_leg=extraleg)
    writefig(filename ,folder='output', write_pdf=False, write_png=True)
    assert False


    ##########################################
    ### Plot simulation against experiment ###
    ##########################################
    
    # load experimental reference - save roughness values (capacitance & curren based) + potential vs. selectivity
    rgh, d_cupd = load_Ac_alloy_data("Ji_COR_CuPd", -1.4)
    for k in rgh:
        print(k, rgh[k], rgh[k][3]*0.4) #rgh[k][2]*2.5)
    okeys = ['CuPd', 'd-CuPd', 'Cu3.4Pd', 'Cu0.3Pd', 'Cu-NP'] 
    dat = {kf(k):d_cupd[k] for k in okeys if k != 'Cu-NP'}
    # load CuAg data
    #akeys = ['CuAg_0%', 'CuAg_1%']
    akeys = ['CuAg_0%', 'CuAg_0.25%', 'CuAg_0.5%', 'CuAg_1%'][::-1]
    r_cuag, d_cuag = load_Ac_alloy_data("Sargent_COR_CuAg", -1.6, ulim=300)
    for k in r_cuag:
        print(k, r_cuag[k], r_cuag[k][3]*5)
    dat.update({kf(k):d_cuag[k] for k in akeys if k != 'CuAg_0%'})
    dat.update({kf('Cu-NP'):d_cupd['Cu-NP']}); dat.update({kf('CuAg_0%'):d_cuag['CuAg_0%']})
    
    cls = plt.cm.cool(np.linspace(0,1,4))
    ls_args = {kf(okeys[i]):dict(ls='--', marker='x', markersize=7, markeredgewidth=2, color=cls[i]) for i in range(len(okeys)-1)}

    cls2 = plt.cm.summer(np.linspace(0,1,4))
    #cls2 = plt.cm.Wistia(np.linspace(0,1,3))
    ls_args.update({kf(akeys[i]):dict(ls='--', marker='o', color=cls2[i]) for i in range(len(akeys)-1)})
    for i in [1,2]:
        ls_args[kf(akeys[i])].update({'alpha':0.2})
    ls_args.update({'Cu-NP':dict(ls='--', marker='x', markersize=7, markeredgewidth=2, color='orange')})
    ls_args.update({r'Cu$_2$O':dict(ls='--', marker='o', color='tab:red')})
  
    # to NOT inclue simulated data in explicitly plotted data set
  # dat.update(sim_dat)
  # kss = list(sim_dat.keys())
  # ls_args.update({kss[i]:dict(ls='-', color=cls[i]) for i in range(len(kss))})
    
    for k in dat:
        dat[k][:,1] *= 100.
    for k in sim_dat:
        sim_dat[k][:,1] *= 100.


    xlabel = r'$\rm U_{SHE}~(V)$'; ylabel = r'selectivity Ac / C$_2$ (%)'
  # plot_xy("ref_CORR_CuPd_SelAc", dat, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, title=r'CORR on Cu/CuPd')


    filename = "ref_CORR_CuPd_SelAc"
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.16,0.98,0.17,0.92],fsize=10)
    ax = plot_xy_ax(deepcopy(dat), ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=6, \
        legpos=1, title=r'')#, extra_leg=extraleg)
        #legpos=1, title=r'CORR on Cu/CuPd')#, extra_leg=extraleg)
    #ax.annotate("", xy=(2.5, 10), xytext=(2.5, 26), arrowprops={"arrowstyle":'<->', "color":'darkgray'})
    #ax.fill_between(x=Us, y1=sim_dat[r'$\rho = 10$'][:,1], y2=sim_dat[r'$\rho = 35$'][:,1], \
    ax.fill_between(x=Us, y1=sim_dat[r'$\rho = 10$'][:,1], y2=sim_dat[r'$\rho = 50$'][:,1], \
        color=clrs['lightblue'], alpha=0.3, zorder=0)
    ax.annotate(r"$\rho=10-50$", xy=(-1.78, 87), size=7, color=clrs['darkblue'])
    ax.fill_between(x=Us, y1=sim_dat[r'$\rho = 100$'][:,1], y2=sim_dat[r'$\rho = 200$'][:,1], \
        color=clrs['orange'], alpha=0.3, zorder=0)
    ax.annotate(r"$\rho=100-200$", xy=(-1.78, 15), size=7, color=clrs['darkred'])
    ax.set_xlim(-1.83, -1.205)
    #ax.set_xlim(-1.83, -1.22)
    ax.set_ylim(1.5974113036320725, 97.91947845208071) #hardcoded
    ax.annotate('(a)', xy=(0.04, 0.95), xycoords='figure fraction')
    writefig(filename ,folder='output', write_pdf=False, write_png=True)

  # assert False
    
    ######################
    ### roughness plot ###
    ######################

    # experimental data
    # adjust relative rgh_i
    for k in rgh:
        rgh[k][3] *= 0.4
    for k in r_cuag:
        r_cuag[k][3] *= 5.0
    rgh.update(r_cuag)

    #for v in [-1.7, -1.6, -1.5]:
    v = -1.7
    y_pot = cluster_pot(dat, v)
    
    rgh_c = {kf(k):rgh[k][0] for k in rgh if k not in akeys}
    #rgh_i = {kf(k):rgh[k][2]*2.5 for k in rgh}
    rgh_i = {kf(k):rgh[k][3] for k in rgh} # rel multiplicator befor merge rgh

    dat_c = {k:np.array([[rgh_c[k], y_pot[k]]]) for k in rgh_c}
    dat_i = {k:np.array([[rgh_i[k], y_pot[k]]]) for k in rgh_i}
    
    # simulation data
    rghs = [1,2,3,4,5,8,10] + np.arange(10,300,10).tolist()
    Us = np.array([v])
    dat_sim = potsweep_rgh_pH(Us=Us, pHs=pHs, \
        rghs=rghs, Lx=[0.15e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pol_acetate_gpaw', tol_diff=1e-8)
    dat_sim = np.array([[float(k),dat_sim[k][pHs[0]]['SAcC2'][0]*100.] for k in dat_sim])


    # plotting
    plot_SI = True
    #plot_SI = False

    filename = "ref_CORR_CuPd_SelAc_rgh_%.1f"%v
 ###filename = "ref_CORR_CuPd_SelAc_rgh_%.1f_a"%v
 ###for k in ['Cu$_2$O', 'CuAg$_{0.25\\%}$', 'CuAg$_{0.5\\%}$', 'CuAg$_{1\\%}$']:
 ###    del dat_i[k]

    if plot_SI:
        filename += "_SI"
    xlabel = r'roughness $\rho$'
    for k in ls_args:
        ls_args[k]['ls'] = 'None' 
    ls_args['Cu-NP']['color'] = 'orange'
    for i in [1,2]:
        del ls_args[kf(akeys[i])]['alpha']
    #extraleg = {k:dict(marker='s', ls='None', color=ls_args[k]['color']) for k in ls_args}
    extraleg = {k:ls_args[k] for k in ls_args}
    extraleg.update({
       #r'$\rho_C^{\mathrm{foil}}$':dict(ls='None', marker='o', color='k'),\
       #r'$\rho_I^{\mathrm{min}}$':dict(ls='None', marker='x', color='k'),\
        'model':dict(ls='-', color='r')})
    
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.16,0.98,0.18,0.92],fsize=10)
    ax = plot_xy_ax(dat_i, ylabel, xlabel, ls_args, tag='', line=[], fsize=10, lsize=6, \
        legpos=1, legcol=3, extra_leg=extraleg)
        #legpos=1, legcol=3, title=r'CORR on Cu/CuPd @ %.1f V$_{\mathrm{SHE}}$'%v, extra_leg=extraleg)
    if plot_SI:
        for k in dat_c:
            ax.plot(dat_c[k][:,0], dat_c[k][:,1], marker='x', color=ls_args[k]['color'], alpha=0.1)
            ax.annotate('', xy=dat_i[k][0], xytext=dat_c[k][0], arrowprops={"arrowstyle":'->', "color":ls_args[k]['color'], 'alpha':0.1})
    anargs = dict(xytext=(380,35), alpha=0.2, size=7, arrowprops={"arrowstyle":'->', "color":'k', 'alpha':0.2})
    if plot_SI:
        ax.annotate(r'$\rho_C^{\mathrm{foil}}$', xy=(380,65), **anargs)
        ax.annotate(r'$\rho_C^{\mathrm{foil}}$', xy=(310,25), **anargs)
    
    anargs = dict(xytext=(2,25), alpha=0.5, size=7, arrowprops={"arrowstyle":'->', "color":'k', 'alpha':0.5})
    if plot_SI:
        ax.annotate(r'$\rho_I^{\mathrm{rel}}$', xy=(15,45), **anargs)
        ax.annotate(r'$\rho_I^{\mathrm{rel}}$', xy=(160,22), **anargs)
        ax.annotate('(b)', xy=(0.04, 0.95), xycoords='figure fraction')
    
    if not plot_SI:
        ax.set_xlim(ax.get_xlim()[0], 310)
        t = ax.text(130, 40, "alloying",
            ha="center", va="center", rotation=-25, size=8,
            bbox=dict(boxstyle="larrow,pad=0.3",
                      fc="none", ec="k", lw=1))
        add_sketch(ax, r'CO', r'C$_{2}$', r'  Ac', dxy=(0.2, -0.1))
    #ax.annotate(r'@ %.1f V$_{\mathrm{SHE}}$'%v, xy=(0.8,1.03), xycoords='axes fraction', size=7)
    
    ax.plot(dat_sim[:,0], dat_sim[:,1], color='r')
    ax.set_ylim(11.630705628121804, 100.36820853571174) # for VJB buil-up

    writefig(filename ,folder='output', write_pdf=False, write_png=True)

