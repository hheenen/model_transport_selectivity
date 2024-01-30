#!/usr/bin/env python

from header import *


def load_literature_CO2RR_COsel_data():
    """
      helper-function to load Pt-ORR data

    """
    fkeys= ["Huang-111", "Huang-OD-Cu", "Kanan-pc-Cu", "Kuhl-pc-Cu"]
    selCO = {fkey:np.loadtxt("../literature_data/dat_CO2RR_CO_%s.txt"%fkey)[:,[0,-1]]
        for fkey in fkeys}
    for k in selCO:
        selCO[k][:,1] *= 100.
    
    # rgh vals as described in SI -- decimals are to identify coloring
    rghvals = {'Huang-OD-Cu':32, 'Kuhl-pc-Cu':1.1,
        'Huang-111':1.2, 'Kanan-pc-Cu':30.1}

    return selCO, rghvals


def plot_CO2RR_CO_pot(filename, dat, ls_args, extraleg):
    """
      helper-function to plot ORR data

    """
    set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05, \
        lrbt=[0.17,0.98,0.18,0.92],fsize=10)
    
    xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity CO / C$_{1} \! + \! \mathrm{C}_{2\!+} \! (\%)$'
    
    ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], lsize=7, \
        legpos=4, title='', extra_leg=extraleg)

    ax.annotate('(a)', xy=(0.08, 0.95), xycoords='figure fraction')
    clr1 = ls_args['sim_1']['color']
    clr30 = ls_args['sim_30']['color']
    ax.annotate(r'$\rho$=1', xytext=(-1.62, 25), xy=(-1.45,25), color=clr1, size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clr1))
    ax.annotate(r'$\rho$=30', xytext=(-1.25, 5), xy=(-1.35,10), color=clr30, size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clr30))
    
    add_sketch2(ax, r'CO$_2$', r'C$_{1} \! + \! \mathrm{C}_{2\!+}$', r'  CO') #dxy=(0.22,-0.2))

    ax.set_xlim(-1.70252, -0.65) #hardcoded
    ax.set_ylim(-4.93, 104.99) #hardcoded
    
    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=False)


def make_plot_CO2RR_CO_pot():
    """
      helper-function to make full plot

    """
    # run model for CO in CO2RR on Cu example
    out_plot = {}
    for rgh in [1,30]:
        sdict['rghs'] = [rgh]
        out_plot.update(run_model_for_example({'sim_%i'%rgh:cdict}, \
            {'sim_%i'%rgh:sdict}, okey='U_SHE'))
        out_plot['sim_%i'%rgh][:,1] *= 100.
        # remove non-converged points
        ind = (np.isnan(out_plot['sim_%i'%rgh][:,1]) == False)
        out_plot['sim_%i'%rgh] = out_plot['sim_%i'%rgh][ind,:]
    
    # load literature data
    selCO, rghCO = load_literature_CO2RR_COsel_data()
    
    # prepare ls_args 
    clsk = {1:clrs['darkblue'], 1.1:clrs['azurblue'], 1.2:clrs['lightblue2'], 
        1.3:clrs['lightblue3'], 1.4:'b',
        30:'r', 30.1:clrs['orange'], 32:clrs['darkyellow'],
        198:clrs['darkred'], 475.0:clrs['darkred']}

    ls_args = {k:dict(ls='--', marker='x', color=clsk[rghCO[k]]) for k in selCO}
    for rgh in [1,30]:
        ls_args.update({'sim_%i'%rgh:dict(ls='-', lw=2.5, color=clsk[rgh])})

    # prepare extra-legend
    extraleg = {r'Cu(111) $\rho=$1':dict(ls='--', marker='x', color=clsk[1.2]),
        r'pc-Cu $\rho=$1':dict(ls='--', marker='x', color=clsk[1.1]),
        r'OD-Cu $\rho=$30':dict(ls='--', marker='x', color=clsk[30.1]),
        r'OD-Cu $\rho=$32':dict(ls='--', marker='x', color=clsk[32]),
        r'model':dict(ls='-', color='r'),
    }

    # make plot
    filename = "Fig3a_CO2RR_Cu_SelCO_pot"
    selCO.update(out_plot)
    plot_CO2RR_CO_pot(filename, selCO, ls_args, extraleg)


#################################
### Parameters for simulation ###
#################################

# dG --> G_des - G_ads (i.e. positive means adsorption is exothermic)
ddG_des = 0.165     # guess

# parameters for CO on Cu
cH_CO = 1100        # Henry's constant for CO atm/M  (10.5194/acp-15-4399-2015)
cD_CO = 20.3e-10    # m2/s              (Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems)
dG_CO_111 = 0.1     # (10.1039/C0EE00071J
dG_red_CO = 1.5     # changed for CO selectivity to 1.5

# effectively used model parameters
cdict = {'cD':[cD_CO], 'cH':[cH_CO], 'dG':[dG_CO_111], \
    'ddG_des':[ddG_des], 'ddG_red':[dG_red_CO]}

# calculation instructions
sdict = {"Us":np.arange(-1.6,-0.6,0.05), "rghs":None, "mdls":[1]}


if __name__ == "__main__":
    make_plot_CO2RR_CO_pot()


