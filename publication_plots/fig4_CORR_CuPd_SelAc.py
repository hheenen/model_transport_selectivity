#!/usr/bin/env python

from header import *
from model_Sel_CORAc.simtools import potsweep_rgh_pH, make_facet_mkin_model_gpaw 
from model_Sel_CORAc.transport.flux_conversion import f2j


def load_literature_CO2RR_COsel_data():
    """
      helper-function to load Pt-ORR data

    """
    tag_dict = {'Cu-Flower':r"Cu-flower $\rho$=390",
        'pc-Cu':r"pc-Cu $\rho$=1",
        'OD-Cu':r"OD-Cu $\rho$=87",
    }
    fkeys = ["d-CuPd", "CuPd", "Cu0.3Pd", "Cu3.4Pd", \
             "CuAg_1%", "CuAg_0.5%", "CuAg_0.25%", \
             "Cu-NP", "CuAg_0%"]
    selAc = {kf(fkey):np.loadtxt(\
        "../literature_data/dat_CO2RR_CO_%s.txt"%fkey)[:,[0,-2]] \
        for fkey in fkeys}
    errAc = {kf(fkey):np.loadtxt(\
        "../literature_data/dat_CO2RR_CO_%s.txt"%fkey)[:,[0,-1]] \
        for fkey in fkeys}
    for k in selAc:
        selAc[k][:,1] *= 100.
        errAc[k][:,1] *= 100.
        errAc[k][:,0] = 0. # x-error not taken
        # unreasonable errors from error-propagation > 50% are ignored:
        ind50 = np.where(np.absolute(errAc[k][:,1]/selAc[k][:,1]) > 0.5)[0]
        errAc[k][ind50,1] = 0.
   
    return selAc, errAc


def kf(k):
    """
      helper-function to create labels

    """
    if k == 'Cu3.4Pd':
        return(r'CuPd$_{23\%}$')
    elif k == 'Cu0.3Pd':
        return(r'CuPd$_{77\%}$')
    elif k == 'CuPd':
        return(r'o-CuPd')
    elif k == 'd-CuPd':
        return(r'CuPd$_{50\%}$')
    elif k == 'CuAg_0%':
        return r'Cu$_2$O'
    elif k == 'CuAg_0.25%':
        return r'CuAg$_{0.25\%}$'
    elif k == 'CuAg_0.5%':
        return r'CuAg$_{0.5\%}$'
    elif k == 'CuAg_1%':
        return r'CuAg$_{1\%}$'
    else:
        return k
    

def cluster_pot(dat, v):
    """
      helper-function to collect potentials at same potential

    """
    sel_v = {k: dat[k][np.absolute(dat[k][:,0] - v).argmin(),1] for k in dat}
    ind_v = {k: np.absolute(dat[k][:,0] - v).argmin() for k in dat}
    # determine standard deviation from target potential "v"
    dif_v = {k: np.absolute(dat[k][:,0] - v).min() for k in dat}
    std_v = np.sqrt(sum([dif_v[k]**2.0 for k in dif_v]) / len(dif_v))
    print("deviation potential dependent data-points: +/- %.3f V"%std_v)
    return sel_v, ind_v


def run_acetate_model_potential():
    """
      helper-function to run acetate model for Fig S4a

    """
    mkin = make_facet_mkin_model_gpaw(100)
    rghs = [10, 50, 100, 200]
    pHs=[14.0]
    Us = np.arange(-1.9, -1.05, 0.1)

    dat = potsweep_rgh_pH(Us=Us, pHs=pHs, \
        rghs=rghs, Lx=[0.15e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pol_acetate_gpaw', tol_diff=1e-8)
    
    sim_dat = {
        r'$\rho = 10$':np.array([dat[10][14.0]['v'], dat[10][14.0]['SAcC2']]).T,
        r'$\rho = 50$':np.array([dat[50][14.0]['v'], dat[50][14.0]['SAcC2']]).T,
        r'$\rho = 100$':np.array([dat[100][14.0]['v'], dat[100][14.0]['SAcC2']]).T,
        r'$\rho = 200$':np.array([dat[200][14.0]['v'], dat[200][14.0]['SAcC2']]).T,
    }
    
    for k in sim_dat:
        sim_dat[k][:,1] *= 100.
    
    return sim_dat


def plot_CO2RR_Acdh_pot(filename, dat, sim_dat, ls_args, derr={}):
    """
      helper-function to plot acetate in CORR on Cu/Pd data

    """
    set_plotting_env(width=3.37,height=3.37/golden_ratio *1.05, \
        lrbt=[0.16,0.98,0.17,0.92],fsize=10)

    xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity Ac / C$_{2\!+} \! (\%)$'

    Us = sim_dat[list(sim_dat.keys())[0]][:,0]

    # plot data-points
    ax = plot_xy_ax(deepcopy(dat), ylabel, xlabel, ls_args, tag="", line=[], lsize=6, \
        legpos=1, title=r'', derr=derr)

    # plot sim-dat as ranges
    ax.fill_between(x=Us, y1=sim_dat[r'$\rho = 10$'][:,1], y2=sim_dat[r'$\rho = 50$'][:,1], \
        color=clrs['lightblue'], alpha=0.3, zorder=0)
    ax.annotate(r"$\rho=10-50$", xy=(-1.78, 87), size=7, color=clrs['darkblue'])
    ax.fill_between(x=Us, y1=sim_dat[r'$\rho = 100$'][:,1], y2=sim_dat[r'$\rho = 200$'][:,1], \
        color=clrs['orange'], alpha=0.3, zorder=0)
    ax.annotate(r"$\rho=100-200$", xy=(-1.78, 15), size=7, color=clrs['darkred'])
    
    # format figure
    ax.set_xlim(-1.83, -1.205)
    #ax.set_xlim(-1.83, -1.22)
    ax.set_ylim(1.5974113036320725, 97.91947845208071) #hardcoded
    ax.annotate('(a)', xy=(0.04, 0.95), xycoords='figure fraction')
    
    writefig(filename ,folder='output', write_pdf=False, write_png=True)


def make_plot_CO2RR_Ac_pot(no_errorbars=True):
    """
      helper-function to make full plot S4a

    """
    # run model for Ac in CORR on Cu/Pd example - potential range
    sim_dat = run_acetate_model_potential()

    # load experimental reference data
    dat, err = load_literature_CO2RR_COsel_data()
    if no_errorbars: # only plot on choice for visibility
        err = {}

    # prepare ls_args
    ls_args = make_lsargs()

    filename = "FigS4a_CORR_CuPd_SelAc_pot"
    plot_CO2RR_Acdh_pot(filename, dat, sim_dat, ls_args, derr=err)


def plot_CO2RR_Acdh_rgh(filename, dat_i, dat_c, sim_dat, ls_args, derr, plot_SI):
    """
      helper-function to plot acetate in CORR on Cu/Pd data

    """
        
    set_plotting_env(width=3.25,height=3.25/golden_ratio *1.05, \
        lrbt=[0.16,0.98,0.185,0.92],fsize=10)
    
    # make extra-legend
    extraleg = {k:ls_args[k] for k in ls_args}
    extraleg.update({
        'model':dict(ls='-', color='r')})
    extraleg.update({' ':dict(ls=None, color='w')})
    extraleg.update({'  ':dict(ls=None, color='w')})
    korder = ['Cu-NP', 'Cu$_2$O', 'model', ' ', 'CuAg$_{1\\%}$', \
        'CuAg$_{0.5\\%}$', 'CuAg$_{0.25\\%}$', '  ', 'o-CuPd', 'CuPd$_{50\\%}$', \
        'CuPd$_{23\\%}$', 'CuPd$_{77\\%}$'] 
    extraleg = {k:extraleg[k] for k in korder}
    
    xlabel = r'roughness $\rho$'; ylabel = r'selectivity Ac / C$_{2\!+} \! (\%)$'
    
    # actual plotting
    ax = plot_xy_ax(dat_i, ylabel, xlabel, ls_args, tag='', line=[], lsize=6, \
        legpos=1, legcol=3, csp=0.8, extra_leg=extraleg, derr=derr)
    
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
        t = ax.text(105, 25, "alloying",
            ha="center", va="center", rotation=0, size=8,
            bbox=dict(boxstyle="larrow,pad=0.3",
                      fc="none", ec="k", lw=1))
        add_sketch(ax, r'CO', r'C$_{2\!+}$', r'  Ac', dxy=(0.2, -0.2))
    
    ax.plot(sim_dat[:,0], sim_dat[:,1], color='r')
    ax.set_ylim(11.630705628121804, 100.36820853571174) # for VJB buil-up
 
    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=False)


def make_lsargs():
    """
      helper-function to make full plot S4a

    """
    okeys = ['CuPd', 'd-CuPd', 'Cu3.4Pd', 'Cu0.3Pd', 'Cu-NP'] 
    cls = plt.cm.cool(np.linspace(0,1,4))
    ls_args = {kf(okeys[i]):dict(ls='--', marker='x', markersize=7, markeredgewidth=2, color=cls[i]) for i in range(len(okeys)-1)}

    akeys = ['CuAg_0%', 'CuAg_0.25%', 'CuAg_0.5%', 'CuAg_1%'][::-1]
    cls2 = plt.cm.summer(np.linspace(0,1,4))
    ls_args.update({kf(akeys[i]):dict(ls='--', marker='o', color=cls2[i]) for i in range(len(akeys)-1)})
    
    for i in [1,2]:
        ls_args[kf(akeys[i])].update({'alpha':0.2})
    ls_args.update({'Cu-NP':dict(ls='--', marker='x', markersize=7, markeredgewidth=2, color='orange')})
    ls_args.update({r'Cu$_2$O':dict(ls='--', marker='o', color='tab:red')})

    return ls_args


def make_plot_CO2RR_Ac_rgh(no_errorbars=True):
    """
      helper-function to make full plot 4 and S4b

    """
    # load experimental reference data
    dat, err = load_literature_CO2RR_COsel_data()
    if no_errorbars: # only plot on choice for visibility
        err = {}
    
    rgh_i = {'Cu-NP': 157.22, 'CuPd': 43.82, 'd-CuPd': 27.31, 'Cu3.4Pd': 75.35, 
        'Cu0.3Pd': 1.20, 'CuAg_0%': 188.20, 'CuAg_0.25%': 19.12, 
        'CuAg_0.5%': 60.97, 'CuAg_1%': 5.0}
    rgh_c = {'Cu-NP': 305.1, 'CuPd': 89.3, 'd-CuPd': 353.1, 'Cu3.4Pd': 357.9, 'Cu0.3Pd': 419.3}

    rgh_i = {kf(k):rgh_i[k] for k in rgh_i}
    rgh_c = {kf(k):rgh_c[k] for k in rgh_c}

    v = -1.7 # potential at which to plot
    y_pot, ind = cluster_pot(dat, v)
    y_err = {k:np.array([[0.0, err[k][ind[k],1]]]) for k in err}
    
    dat_c = {k:np.array([[rgh_c[k], y_pot[k]]]) for k in rgh_c}
    dat_i = {k:np.array([[rgh_i[k], y_pot[k]]]) for k in rgh_i}
    
    # simulation data
    mkin = make_facet_mkin_model_gpaw(100)
    rghs = [1,2,3,4,5,8,10] + np.arange(10,300,10).tolist()
    Us = np.array([v])
    sim_dat = potsweep_rgh_pH(Us=Us, pHs=[14.0], \
        rghs=rghs, Lx=[0.15e-4, 1e-4], diffmode='diff-rct-num', mkin=mkin, savekey='pol_acetate_gpaw', tol_diff=1e-8)
    sim_dat = np.array([[float(k),sim_dat[k][14.0]['SAcC2'][0]*100.] for k in sim_dat])

    # prepare ls_args
    ls_args = make_lsargs()
    for k in ls_args:
        ls_args[k]['ls'] = 'None' 
    ls_args['Cu-NP']['color'] = 'orange'
   #for i in [1,2]:
   #    del ls_args[kf(akeys[i])]['alpha']
    
    # plot Fig 4
    filename = "Fig4_CORR_CuPd_SelAc_rgh_%.1f"%v
    plot_CO2RR_Acdh_rgh(filename, dat_i, dat_c, sim_dat, ls_args, derr=y_err, plot_SI=False)

    
    # plot Fig 4
    filename = "FigS4b_CORR_CuPd_SelAc_rgh_%.1f"%v
    plot_CO2RR_Acdh_rgh(filename, dat_i, dat_c, sim_dat, ls_args, derr=y_err, plot_SI=True)



if __name__ == "__main__":
    
    make_plot_CO2RR_Ac_rgh(no_errorbars=True)

    make_plot_CO2RR_Ac_pot(no_errorbars=True)
    
