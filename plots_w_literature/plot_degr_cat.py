#!/usr/bin/env python

from header import *

import pandas as pd 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

if __name__ == "__main__":


    #############
    ### Cu-NP ###
    #############

    
    ###########################################################
    ### read data and transform to dictionary as old format ###
    ###########################################################

    df = pd.read_csv('../experimental_reference/buonsanti.csv')
    d_sCO_t = {}; d_rgh_t = {}; d_sCO_rgh = {}
    for NPs in [16, 41, 65]:
        d_sCO_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], df['S_CO %i'%NPs]*100.]).T})
        rgh_r = df['J C2H4 %i'%NPs]/df['J C2H4 16'][0]
        #rgh_r = df['J Total %i'%NPs]/df['J Total 16'][0]
        d_rgh_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], rgh_r]).T})
        d_sCO_rgh.update({"NP %inm"%NPs:np.array([rgh_r, df['S_CO %i'%NPs]*100.]).T})
       #d_rgh_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], df['Rough %i'%NPs]]).T})
       #d_sCO_rgh.update({"NP %inm"%NPs:np.array([df['Rough %i'%NPs], df['S_CO %i'%NPs]]).T})
    

    ##########################
    ### shared x-axis plot ### 
    ##########################

    filename = "exp_CO2RR_Cu_degr"
    xlabel = r'operation time (h)'
    #ylabel = r'selectivity CO / C$_{1\!+\!2}$ (%)'
    ylabel = r'selectivity CO / C$_{1} \! + \! \mathrm{C}_{2\!+} \! (\%)$'# + 2\!\!+}$ (%)'

    clsk = {16:clrs['darkblue'], 65:clrs['orange'], 41:clrs['darkyellow']}
    ls_args = {"NP %inm"%i:dict(ls='--', marker='x', color=clsk[i]) for i in [16,41,65]}

    
    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.5, lrbt=[0.16,0.98,0.12,0.98],fsize=10)
    fig, axes = plt.subplots(2,1)
    axes[0] = plot_xy_ax_in(axes[0], d_sCO_t, ylabel, "", ls_args, tag="", line=[], fsize=10, lsize=7, \
        legpos=2, title='', nolegend=True)
        #legpos=4, title=r'CO$_2$RR on Cu', extra_leg=extraleg)
    axes[1] = plot_xy_ax_in(axes[1], d_rgh_t, r"rel. roughness $\rho$", xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, \
        legpos=1, title='')#, extra_leg=extraleg)

    [ax.set_xlim(-0.5, 12.5) for ax in axes]
    #plt.xticks([0,2,4,6,8,10,12])
    axes[0].set_xticks([0,2,4,6,8,10,12])
    axes[1].set_xticks([0,2,4,6,8,10,12])
    axes[0].set_xticklabels([])
    plt.subplots_adjust(hspace=0.05)
    axes[0].annotate('(a)', xy=(-0.08, 0.95), xycoords='axes fraction')
    axes[1].annotate('(b)', xy=(-0.08, 0.95), xycoords='axes fraction')

    axes[0].yaxis.set_label_coords(-0.12,0.4)
    axes[0].yaxis.set_label_coords(-0.12,0.3)

    writefig(filename ,folder='output', write_pdf=False, write_png=True)
    

    ###################
    ### single plot ###
    ###################

    filename = "exp_CO2RR_Cu_degr_rgh_sel"
    xlabel = r'relative roughness $\rho$'
    #ylabel = r'selectivity CO / C$_{1\!+\!2}$ (%)'
    ylabel = r'selectivity CO / C$_{1} \! + \! \mathrm{C}_{2\!+} \! (\%)$'# + 2\!\!+}$ (%)'

    set_plotting_env(width=3.37,height=3.37/ golden_ratio *1.05, lrbt=[0.17,0.98,0.18,0.92],fsize=10)
    ax = plot_xy_ax(d_sCO_rgh, ylabel, xlabel, ls_args, tag="", line=[], fsize=10, lsize=7, \
        legpos=1, title='')#, extra_leg=extraleg)

    ax.annotate("", xy=(0.4, 0.8), xycoords='axes fraction',
            xytext=(0.4, 0.55), size=10,
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=1.4", color=clsk[16], fc='w'), zorder=0)
    ax.annotate("", xy=(0.5, 0.4), xycoords='axes fraction',
            xytext=(0.5, 0.2), size=10,
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=2.0", color=clsk[65], fc='w'), zorder=0)
    ax.annotate("", xy=(0.9, 0.2), xycoords='axes fraction',
            xytext=(0.9, 0.1), size=10,
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=2.0", color=clsk[41], fc='w'), zorder=0)
    
    ax.set_ylim(2, ax.get_ylim()[1])
    ax.set_xlim(ax.get_xlim()[0], 2.7)

    writefig(filename ,folder='output', write_pdf=False, write_png=True)


    assert False

    # NOTE looks like we cannot make a neat plot like the other for the ORR

    ##############
    ### Pd-ORR ###
    ##############
    
    # todo:
    # ( ) read-data
    # ( ) make xy-plots
    # ( ) add shared axes

    
    ###########################################################
    ### read data and transform to dictionary as old format ###
    ###########################################################

    df = pd.read_csv('../experimental_reference/roughness_cycles_Pd.csv')
    print(df.index)
    print(df.columns)
    #print(df)

  # d_sCO_t = {}; d_rgh_t = {}; d_sCO_rgh = {}
  # for NPs in [16, 41, 65]:
  #     d_sCO_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], df['S_CO %i'%NPs]*100.]).T})
  #     rgh_r = df['J C2H4 %i'%NPs]/df['J C2H4 16'][0]
  #     d_rgh_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], rgh_r]).T})
  #     d_sCO_rgh.update({"NP %inm"%NPs:np.array([rgh_r, df['S_CO %i'%NPs]*100.]).T})
  #    #d_rgh_t.update({"NP %inm"%NPs:np.array([df['Time (Hours)'], df['Rough %i'%NPs]]).T})
  #    #d_sCO_rgh.update({"NP %inm"%NPs:np.array([df['Rough %i'%NPs], df['S_CO %i'%NPs]]).T})


  # print(d_sCO_t)
  # print(d_rgh_t)

    assert False

    print(df.index)
    print(df.columns)
    print(df)
    

    df['roughness_relative_smallest 16'] = df['J C2H4 16']/df['J C2H4 16'][0]
    df['roughness_relative_smallest 41'] = df['J C2H4 41']/df['J C2H4 16'][0]
    df['roughness_relative_smallest 65'] = df['J C2H4 65']/df['J C2H4 16'][0]
