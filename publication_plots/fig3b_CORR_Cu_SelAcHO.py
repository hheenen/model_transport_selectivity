#!/usr/bin/env python

from header import *


def pltargs(k):
    """
      helper-function to create ls_args

    """
    clsk = {"1":clrs['darkblue'], "5":clrs['lightblue'], \
        "87":clrs['azurblue'], "390":clrs['lightblue']}
    if k.find('model') != -1:
        ls='-'; marker='None'; zorder=1
    else:
        ls='--'; marker='x'; zorder=99
    for t in clsk:
        if k.find(t) != -1:
            clr = clsk[t]
    return dict(ls=ls, marker=marker, color=clr, zorder=zorder)


def load_literature_CO2RR_COsel_data():
    """
      helper-function to load Pt-ORR data

    """
    tag_dict = {'Cu-Flower':r"Cu-flower $\rho$=390",
        'pc-Cu':r"pc-Cu $\rho$=1",
        'OD-Cu':r"OD-Cu $\rho$=87",
    }
    fkeys= ["Cu-Flower" , "OD-Cu", "pc-Cu"]
    selAcdh = {tag_dict[fkey]:np.loadtxt(\
        "../literature_data/dat_CO2RR_CO_COR@%s.txt"%fkey)[:,[0,-1]] \
        for fkey in fkeys}
    for k in selAcdh:
        selAcdh[k][:,1] *= 100.
   
    return selAcdh


def plot_CO2RR_Acdh_pot(filename, dat, ls_args):
    """
      helper-function to plot Acdh in CORR on Cu data

    """
    set_plotting_env(width=3.25, height=3.25/golden_ratio *1.05, \
        lrbt=[0.16,0.98,0.18,0.92], fsize=10)

    xlabel = r'$U\rm _{SHE}~(V)$'; ylabel = r'selectivity MeCHO / C$_{2\!+} \!(\%)$'

    # prep extra legend
    extraleg = {k:ls_args[k] for k in ls_args if k.find('model') == -1}
    extraleg.update({'model':dict(ls='-', color='k')})
    
    ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], lsize=7, \
        legpos=2, title='', extra_leg=extraleg)

    # annotations
    ax.annotate(r'$\rho$=1', xytext=(-1.31, 22), xy=(-1.26,20), color=clrs['darkblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['darkblue']))
    ax.annotate(r'$\rho$=87', xytext=(-1.23, 4), xy=(-1.15, 2), color=clrs['azurblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['azurblue']))
    ax.annotate(r'$\rho$=390', xytext=(-1.05, -1), xy=(-1.01,8), color=clrs['lightblue'], size=8, 
            arrowprops=dict(arrowstyle="-", connectionstyle="arc3", color=clrs['lightblue']))
    
    add_sketch(ax, r'CO', r'C$_{2\!+}$', r'MeCHO', dxy=(0.12, 0.05))
    
    ax.annotate('(b)', xy=(0.08, 0.95), xycoords='figure fraction')
    ax.set_ylim(-4.87, 62) #hardcoded
    
    writefig(filename ,folder='output', write_pdf=False, \
        write_png=True, write_eps=False) 


def make_plot_CO2RR_Acdh_pot():
    """
      helper-function to make full plot

    """
    
    # run model for Ac in CORR on Cu example
    out_plot = {}
    for rgh in [1,87,390]:
        sdict['rghs'] = [rgh]
        out_plot.update(run_model_for_example({'Acdh_%i'%rgh:cdict}, \
            {'Acdh_%i'%rgh:sdict}, okey='U_SHE'))
        out_plot['Acdh_%i'%rgh][:,1] *= 100.
        # remove non-converged points
        ind = (np.isnan(out_plot['Acdh_%i'%rgh][:,1]) == False)
        out_plot['Acdh_%i'%rgh] = out_plot['Acdh_%i'%rgh][ind,:]

    
    # load literature data
    dat = load_literature_CO2RR_COsel_data()

    # add simulation data
    dat.update({r'model $\rho=1$':out_plot['Acdh_1']})
    dat.update({r'model $\rho=87$':out_plot['Acdh_87']})
    dat.update({r'model $\rho=390$':out_plot['Acdh_390']})
    

    # prepare ls_args 
    ls_args = {}
    for k in dat:
        ls_args.update({k:pltargs(k)})

    # make plot
    filename = "Fig3b_CORR_Cu_SelAcdh_pot"
    plot_CO2RR_Acdh_pot(filename, dat, ls_args)


#################################
### Parameters for simulation ###
#################################

# dG --> G_des - G_ads (i.e. positive means adsorption is exothermic)
ddG_des = 0.165      # guess
ddGred_Acdh = 0.96 

# Acetaldehyde on Cu
cD_Acdh = 13.75e-10 # m2/s (10.3390/atmos11101057 average of two)
cH_Acdh = 0.0759    # atm /mol (10.1016/S0045-6535(00)00505-1)
dG_Acdh = -0.14 

# effectively used model parameters
cdict = {'cD':[cD_Acdh], 'cH':[cH_Acdh], 'dG':[dG_Acdh], \
    'ddG_des':[ddG_des], 'ddG_red':[ddGred_Acdh]}

# calculation instructions
sdict = {"Us":np.arange(-1.35,-1.0,0.05), "rghs":None, "mdls":[1]}
         

if __name__ == "__main__":
    make_plot_CO2RR_Acdh_pot()


