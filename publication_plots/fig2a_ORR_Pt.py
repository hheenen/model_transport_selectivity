#!/usr/bin/env python

from header import *
            

def load_literature_ORR_Pt_data():
    """
      helper-function to load Pt-ORR data

    """
    dat1 = np.loadtxt("../literature_data/dat_ORR_Pt_Pc.txt")
    dat2 = np.loadtxt("../literature_data/dat_ORR_Pt_22.txt")
    dat3 = np.loadtxt("../literature_data/dat_ORR_Pt_10.txt")
    rgh = [1., 0.22, 0.1]

    # match potentials of different measurements and take average voltage
    # - so rho can be compared; max std = 0.012 V_RHE
    ind2 = []; ind3 = []
    for v in dat1[:,0]:
        ind2.append(np.absolute(dat2[:,0] - v).argmin())
        ind3.append(np.absolute(dat3[:,0] - v).argmin())
    dat2 = dat2[ind2,:]; dat3 = dat3[ind3,:]
    v = np.array([dat1[:,0], dat2[:,0], dat3[:,0]]).T
    v = np.around(np.mean(np.array([dat1[:,0], dat2[:,0], dat3[:,0]]).T, axis=1),4)
    #print("error", np.std(np.array([dat1[:,0], dat2[:,0], dat3[:,0]]).T, axis=1).max())

    # save selected potentials - otherwise too many
    out = {}
    ind = [8, 9, 10, 11]
    tag = r'%.2f V vs. RHE'
    for ii in ind:
        out.update({tag%v[ii]:np.array([rgh, [dat1[ii,1], dat2[ii,1], dat3[ii,1]]]).T})
    
    return out


def plot_ORR_Pt(filename, dat, ls_args, sketch):
    """
      helper-function to plot ORR data

    """
        
    set_plotting_env(width=3.25,height=3.25/ golden_ratio *1.05,\
               lrbt=[0.16,0.98,0.185,0.92],fsize=10)
    
    xlabel = r'roughness $\rho$'; ylabel = r'selectivity H$_2$O$_2$ (%)'
    
    ax = plot_xy_ax(dat, ylabel, xlabel, ls_args, tag="", line=[], lsize=7)#, title=r'ORR on Pt')
    ax.set_xlim(0.055, 1.045) # hardcoded if model is not included
    ax.set_ylim(3.7159276980511073, 35.53675756282097) # hardcoded -- no model
    
    # both sub-plots are (a)
    ax.annotate('(a)', xy=(0.04, 0.95), xycoords='figure fraction')

    if sketch:
        add_sketch(ax, r'O$_2$', r'H$_2$O', r'H$_2$O$_2$')
 
    writefig(filename ,folder='output', write_pdf=False, write_png=True, write_eps=False) 


def make_plots_ORR_Pt():
    """
      helper-function to make full plots

    """

    # run model for ORR on Pt example
    out_plot = run_model_for_example({'H2O2_Pt':cdict}, {'H2O2_Pt':sdict})
    # out_plot['H2O2_Pt'][:,0] /= out_plot['Fmdh'][:,0].min() # relative roughness
    out_plot['H2O2_Pt'][:,1] *= 100.                        # convert percentage
    
    # load literature data
    selORR = load_literature_ORR_Pt_data()
    
    # prepare plotting colors
    ks = list(selORR.keys())
    cls = plt.cm.jet(np.linspace(0,1,len(ks)))
    #ls_args = {r'ORR@Cu $^{[]}$':dict(ls='--', marker='x', color='k'), 
    ls_args = {ks[i]:dict(ls='--', marker='x', color=cls[i], alpha=0.3) \
        for i in range(len(ks))}
    ls_args[r'0.58 V vs. RHE']['alpha'] = 1.0
    ls_args.update({r'model':dict(ls='-', color='r')})
    ls_args.update({r'pc-Pt':dict(ls='--', marker='x', color='k')})

    # plot Fig2a
    filename = "Fig2a_ORR_Pt_SelH2O2"
    # save r'0.58 V vs. RHE' as pc-Pt for main text
    dat = {r'pc-Pt':selORR[r'0.58 V vs. RHE'], r'model':out_plot['H2O2_Pt']}
    plot_ORR_Pt(filename, dat, ls_args, sketch=True)

    # plot FigS4a
    filename = "FigS4a_ORR_Pt_SelH2O2"
    dat = selORR; dat.update({r'model':out_plot['H2O2_Pt']})
    plot_ORR_Pt(filename, dat, ls_args, sketch=False)


#################################
### Parameters for simulation ###
#################################

# dG --> G_des - G_ads 
ddG_h2o2_des = 0.2   # guess 
dG_red = 1.1         # guess

# H2O2
cD_h2o2 = 15.6e-10     # m2/s --> 10.1039/C2CP40616K various values electrolyte depended
cH_h2o2 = 1.20e-05     # atm / M  --> 10.1021/jp951168n

# used model parameters
cdict = {'cD':[cD_h2o2], 'cH':[cH_h2o2], 'dG':[0.1], \
    'ddG_des':[ddG_h2o2_des], 'ddG_red':[dG_red]}

# calculation parameters
sdict = {"Us":[-1.8], "rghs":np.arange(0.1,1.0,0.05), "mdls":[1]}


if __name__ == "__main__":
    make_plots_ORR_Pt()
    
